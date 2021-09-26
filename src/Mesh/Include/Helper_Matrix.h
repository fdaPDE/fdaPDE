#ifndef __HELPER_MATRIX_H__
#define __HELPER_MATRIX_H__

//Class storing a matrix of "list"
template<typename T>
class HelperMatrix{
public:
    using value_type = std::vector<T>;
    using container_type = std::vector< value_type >;

    HelperMatrix( const container_type& matr, UInt nrows, UInt ncols): matr_(matr),nrows_(nrows),ncols_(ncols){}

    HelperMatrix(UInt nrows, UInt ncols);

    value_type& operator[](UInt j);
    const value_type& operator[] (UInt j) const;

    value_type& operator()(UInt i , UInt j);
    const value_type& operator() (UInt i, UInt j) const;

    const UInt nrows() const { return nrows_;};
    const UInt ncols() const { return ncols_;};

    //debugging
    template <typename U>
    friend std::ostream& operator<< (std::ostream&, const HelperMatrix<U>&);

private:
    container_type matr_;
    const UInt nrows_;
    const UInt ncols_;


};

template< typename T>
HelperMatrix<T>::HelperMatrix(UInt nrows, UInt ncols):nrows_(nrows),ncols_(ncols){

    //fill the container with empty value_type
    //value type must have a default constructor
    matr_ = container_type(nrows_*ncols_, value_type({}));
}

template< typename T>
typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator[](UInt j){
    return matr_[j];
}

template< typename T>
const typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator[] (UInt j) const{
    return matr_[j];
}
template<typename T>
typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator()(UInt i , UInt j){
    return matr_[i+nrows_*j];
}

template<typename T>
const typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator() (UInt i, UInt j) const{
    return matr_[i+nrows_*j];
}

template <typename T>
std::ostream& operator<< (std::ostream& stream_, const HelperMatrix<T>& matr_){
    for(UInt i = 0; i < matr_.nrows_; ++i){
        for(UInt j =0; j< matr_.ncols_; ++j ) {
            stream_ << "(" << i << ", " << j << ")" << " :";
            for (const auto elem : matr_(i, j))
                stream_ << elem << " ";
            stream_<<std::endl;
        }
    }
    return stream_;
}

//RIntegerMatrix specialization
template<>
class HelperMatrix<RIntegerMatrix>{
    public:
        using value_type = RIntegerMatrix;
        using container_type = std::vector< value_type >;

        // Robject must be filled!
        HelperMatrix(SEXP Robject);

        // Robject(i,j) is allocated according to lengths
        HelperMatrix(SEXP Robject, const RIntegerMatrix& lenghts);

        value_type& operator[](UInt j);
        const value_type& operator[] (UInt j) const;

        value_type& operator()(UInt i , UInt j);
        const value_type& operator() (UInt i, UInt j) const;

        const UInt nrows() const { return nrows_;};
        const UInt ncols() const { return ncols_;};

    private:
        container_type matr_;
        const UInt nrows_;
        const UInt ncols_;
};

HelperMatrix<RIntegerMatrix>::HelperMatrix(SEXP Robject):nrows_( INTEGER(Rf_getAttrib(Robject, R_DimSymbol))[0] ),
                                                         ncols_( INTEGER(Rf_getAttrib(Robject, R_DimSymbol))[1] ){
    matr_.reserve(nrows_*ncols_);
    for(UInt i=0; i<nrows_*ncols_; ++i){
            matr_.emplace_back(VECTOR_ELT(Robject,i));
        }
}

HelperMatrix<RIntegerMatrix>::HelperMatrix(SEXP Robject, const RIntegerMatrix& lengths):nrows_( INTEGER(Rf_getAttrib(Robject, R_DimSymbol))[0] ),
                                                                                        ncols_( INTEGER(Rf_getAttrib(Robject, R_DimSymbol))[1] ){
    matr_.reserve(nrows_*ncols_);
    for(UInt i=0; i<nrows_*ncols_; ++i){
        SET_VECTOR_ELT(Robject,i,Rf_allocMatrix(INTSXP,1,lengths[i]));
        matr_.emplace_back(VECTOR_ELT(Robject,i));
    }
}

HelperMatrix<RIntegerMatrix>::value_type& HelperMatrix<RIntegerMatrix>::operator[](UInt j){
    return matr_[j];
}

const HelperMatrix<RIntegerMatrix>::value_type& HelperMatrix<RIntegerMatrix>::operator[] (UInt j) const{
    return matr_[j];
}

HelperMatrix<RIntegerMatrix>::value_type& HelperMatrix<RIntegerMatrix>::operator()(UInt i , UInt j){
    return matr_[i+nrows_*j];
}

const HelperMatrix<RIntegerMatrix>::value_type& HelperMatrix<RIntegerMatrix>::operator()(UInt i , UInt j)const{
    return matr_[i+nrows_*j];
}

#endif //__HELPER_MATRIX_H__