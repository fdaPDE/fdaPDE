/*
    =======================================================================
    =======================================================================
    ==                                                                   ==
    ==  FSPAI:  Factorized SPAI algorithm to compute a Factorized SParse ==
    ==          Approximate Inverse matrix for symmetric positive        ==
    ==          definite systems.                                        ==
    ==                                                                   ==
    ==  Copyright (C)  2010, 2011 by                                     ==
    ==                 Matous Sedlacek <sedlacek@in.tum.de>              ==
    ==                 Chair of Scientific Computing -- Informatics V    ==
    ==                 Technische Universität München                    ==
    ==                                                                   ==
    ==  This file is part of FSPAI.                                      ==
    ==                                                                   ==
    ==  FSPAI is free software: you can redistribute it and/or           ==
    ==  modify it under the terms of the GNU Lesser General Public       ==
    ==  License as published by the Free Software Foundation, either     ==
    ==  version 3 of the License, or (at your option) any later version. ==
    ==                                                                   ==
    ==  FSPAI is distributed in the hope that it will be useful,         ==
    ==  but WITHOUT ANY WARRANTY; without even the implied warranty of   ==
    ==  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    ==
    ==  GNU Lesser General Public License for more details.              ==
    ==                                                                   ==
    ==  You should have received a copy of the GNU Lesser General        ==
    ==  Public License along with FSPAI.                                 ==
    ==  If not, see <http://www.gnu.org/licenses/>.                      ==
    ==                                                                   ==
    =======================================================================
    =======================================================================
*/

template<class T_Field>
PCG<T_Field>::PCG
(   const ENV_Handler&     env_handler_,
    Comm_Handler<T_Field>* comm_handler_,
    Matrix<T_Field>*       mtx_,
    const double           tol_,
    const unsigned int     maxit_,
    const int              para_max_levels_,
    const double           para_threshold_,
    const double           para_filter_) :
        env_handler(env_handler_),
        comm_handler(comm_handler_),
        mtx(mtx_),
        tol(tol_),
        maxit(maxit_),
        para_max_levels(para_max_levels_),
        para_threshold(para_threshold_),
        para_filter(para_filter_)
{
    x = new T_Field[mtx->my_nbr_cols];
    b = new T_Field[mtx->my_nbr_cols];
}



template<class T_Field>
PCG<T_Field>::~PCG( )
{
    if(x) delete [] x;
    if(b) delete [] b;
}



template<class T_Field> void
PCG<T_Field>::Init_Solution
(   const int dim )
{
    Init_Vec(x, dim);
}



template<class T_Field> void
PCG<T_Field>::Create_RHS
(   const int rhs_choice,
    const int dim )
{
    switch(rhs_choice)
    {
        case 0:  RHS_Rand(); break;
        case 1:  RHS_Ones(); break;
    }
}



template<class T_Field> void
PCG<T_Field>::Set_Precond
(   const Matrix<T_Field>* precond_ )
{
    precond = precond_;
}



template<class T_Field> void
PCG<T_Field>::Solve_System()
{
    // Note that this PCG algorithm is similar to the
    // PCG algorithm from ParaSails. We adapt it to our
    // data structures, modify the communication and
    // make it solve complex systems as well.
    T_Field *p    = NULL,
            *s    = NULL,
            *r    = NULL,
            *ptmp = NULL;
    double   bi_prod,
             i_prod,
             eps,
             gamma,
             gamma_old,
             alpha,
             beta;
    int      dim = mtx->my_nbr_cols;
    unsigned int i = 0;

    Timer timer = Timer();
    // Start time measurement
    timer.Start( env_handler );

    // compute square of absolute stopping threshold
    // bi_prod = <b,b>
    bi_prod = Inner_Prod(dim, b, b );
    eps = (tol*tol)*bi_prod;

    // Check to see if the rhs vector b is zero
    //if (bi_prod == 0.0) //todo: avoid dirty testing against double value
    if (fabs(bi_prod-0.0) < 1.0e-10) // avoiding testing against double value
    {                                // Set x equal to zero and return
        Copy_Vector(dim, b, x);
        return;
    }

    p       = new T_Field[dim]; Init_Vec(p,dim);
    ptmp    = new T_Field[dim]; Init_Vec(ptmp,dim);
    s       = new T_Field[dim]; Init_Vec(s,dim);
    r       = new T_Field[dim]; Init_Vec(r,dim);

    MV_Prod(x, r);                // r = Ax
    Scale_Vector(dim, -1.0, r);   // r = -r
    Axpy(dim, 1.0, b, r);         // r = r + b

    // p = (L_M^T*L_M)*r,
    if (precond != NULL)
    {
        FSPAI_ApplyTrans(ptmp, r);  // ptmp = L_MT*r
        FSPAI_Apply(p, ptmp);       // p = L_M*ptmp;
    }
    else
        Copy_Vector(dim, r, p);

    // gamma = <r,p>
    gamma = Inner_Prod(dim, r, p);

    while( (i+1) <= maxit )
    {
        i++;

        // s = A*p
        MV_Prod(p, s);

        // alpha = gamma / <s,p>
        alpha = gamma / Inner_Prod(dim, s, p);

        gamma_old = gamma;

        // x = x + alpha*p
        Axpy(dim, alpha, p, x);

        // r = r - alpha*s
        Axpy(dim, -alpha, s, r);

        // s = (L_M^T*L_M)*r
        if (precond != NULL)
        {
            FSPAI_ApplyTrans(ptmp, r);  // ptmp = L_MT*r
            FSPAI_Apply(s, ptmp);       // s = L_M*ptmp;
        }
        else
            Copy_Vector(dim, r, s);

        // gamma = <r,s>
        gamma = Inner_Prod(dim, r, s);

        // set i_prod for convergence test
        i_prod = Inner_Prod(dim, r, r);

        // check for convergence
        if (i_prod < eps) break;

        // non-convergence test
        if (i >= 1000 && i_prod/bi_prod > 0.01)
        {
           if (mtx->my_id == 0)
           {
               std::cout << "\n\t    Aborting PCG solve due to" << std::endl;
               std::cout << "\t    slow or no convergence...\t\t ";
           }
           break;
        }
        // beta = gamma / gamma_old
        beta = gamma / gamma_old;

        // p = s + beta p
        Scale_Vector(dim, beta, p);
        Axpy(dim, 1.0, s, p);

     }

     // compute exact relative residual norm
     MV_Prod(x, r);                 // r = Ax
     Scale_Vector(dim, -1.0, r);    // r = -r
     Axpy(dim, 1.0, b, r);          // r = r + b
     i_prod = Inner_Prod(dim, r, r);

     // Stop and report time measurement
     env_handler.Barrier();
     timer.Stop( env_handler );
     timer.Report( env_handler );

     if (mtx->my_id == 0)
     {
         if(i == maxit) std::cout << "\t    Maximum number of iter: "
                                  << i << " reached!" << std::endl;
         else           std::cout << "\t    Converged at iter:\t    "
                                  << i << std::endl;
         std::cout << "\t    With Rel-res-norm:\t    "
                   << sqrt(i_prod/bi_prod) << std::endl;
     }

     delete [] p;
     delete [] s;
     delete [] r;
     delete [] ptmp;
}



template<class T_Field> void
PCG<T_Field>::Solution_To_File
( const char* file )
{
    int     nnz = mtx->n, //dimension of matrix.
            num_procs,
            my_id,
            ierr;
    FILE    *f;
    const char *mm_string = "%%MatrixMarket";
    char    fullname[1024],
            cat_cmd[1024],
            rm_cmd[1024];

    env_handler.Get_Environment_Params( num_procs, my_id );

    if (num_procs > 1)
    {
        sprintf(fullname,   "%s_tmp%5.5d",      file, my_id);
        sprintf(cat_cmd,    "cat %s_tmp* > %s", file, file);
        sprintf(rm_cmd,     "rm -f %s_tmp*",    file);
    }
    else
        sprintf(fullname, "%s", file);

    if ( !( f = fopen( fullname,"w" ) ) )
        throw std::runtime_error(
            "\n\tERROR:  Failed writing solution to file "
            + std::string(fullname) + "\n"
            "\n\t\tCheck your access rights!\n");

    // write Matrix-Market header
    if (my_id == 0)
    {
        fprintf(f, "%s ", mm_string);
        Write_Header(f);
        fprintf(f, "%d 1\n", nnz);
        fflush(f);
    }

    for (int i = 0; i < mtx->my_nbr_cols; i++)
        Write_Line(x[i], f);

    fflush(f);
    fclose(f);

    env_handler.Barrier();

    if (num_procs > 1)
        if (my_id == 0)
            ierr = system(cat_cmd);

    env_handler.Barrier();

    if (num_procs > 1)
        if (my_id == 0)
            ierr = system(rm_cmd);
}

