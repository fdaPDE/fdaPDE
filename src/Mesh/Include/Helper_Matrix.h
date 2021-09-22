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

    UInt size(UInt i )const;
    UInt size(UInt i, UInt j) const;
    UInt size() const;

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

template < typename T>
UInt HelperMatrix<T>::size(UInt i )const{
    return matr_[i].size();
}

template < typename T>
UInt HelperMatrix<T>::size(UInt i, UInt j) const{
    return matr_[i + nrows_*j].size();
}

template< typename T>
UInt HelperMatrix<T>::size() const{
    return nrows_*ncols_;
}

//compute the lengths of each vector in matrix_
template <typename T>
std::vector<UInt> compute_lengths(const HelperMatrix<T>& matrix_){
    std::vector<UInt> lengths(matrix_.size(),0);
    for(UInt i=0; i<matrix_.size(); ++i)
        lengths[i] = matrix_.size(i);
    return lengths;
}

//Fills the neighbors matrix
template <typename T>
void helper_neighbors( HelperMatrix<T>& result, const HelperMatrix<T>& node_list){
    for(UInt node = 0; node < node_list.nrows(); ++node){
        const typename HelperMatrix<T>::value_type& curr_right = node_list(node,0); //questo vector<int> // lista di lati collegati
        // a node da parte "0", i.e da hanno il nodo allora loro dx
        const typename HelperMatrix<T>::value_type& curr_left = node_list(node,1); //questo vector<int> // lista di lati collegati
        // a node da parte "1", i.e da hanno il nodo allora loro sx
        helper_neigh_same_side_(result,curr_right,0);
        helper_neigh_same_side_(result,curr_left,1);

        //filling different sides
        for(const auto& i : curr_right)
            for(const auto& j : curr_left){
                result(i,0).push_back(j);
                result(j,1).push_back(i);
            }
    }
}

//This function fills the neighbors on the same side
template <typename T>
void helper_neigh_same_side_( HelperMatrix<T>& result, const typename HelperMatrix<T>::value_type& curr,UInt idx){

    if( !curr.empty() ) {
        for (UInt i = 0; i < curr.size() - 1; ++i) {
            for (UInt j = i + 1; j < curr.size(); ++j) {
                result(curr[i], idx).push_back(curr[j]);
                result(curr[j], idx).push_back(curr[i]);
            }
        }
    }
}

#endif //__HELPER_MATRIX_H__