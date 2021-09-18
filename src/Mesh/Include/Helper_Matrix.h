#ifndef __HELPER_MATRIX_H__
#define __HELPER_MATRIX_H__

//Class storing a matrix of "list"
template<typename T>
class HelperMatrix{
public:
    using value_type = std::vector<T>;
    using container_type = std::vector< value_type >;

    HelperMatrix( const container_type& matr, int nrows, int ncols): matr_(matr),nrows_(nrows),ncols_(ncols){}

    HelperMatrix(int nrows, int ncols);

    value_type& operator[](int j);
    const value_type& operator[] (int j) const;

    value_type& operator()(int i , int j);
    const value_type& operator() (int i, int j) const;

    const int nrows() const { return nrows_;};
    const int ncols() const { return ncols_;};

    int size(int i )const;
    int size(int i, int j) const;
    int size() const;

    //debugging
    template <typename U>
    friend std::ostream& operator<< (std::ostream&, const HelperMatrix<U>&);

private:
    container_type matr_;
    const int nrows_;
    const int ncols_;


};

template< typename T>
HelperMatrix<T>::HelperMatrix(int nrows, int ncols):nrows_(nrows),ncols_(ncols){

    //fill the container with empty value_type
    //value type must have a default constructor
    matr_ = container_type(nrows_*ncols_, value_type({}));
}

template< typename T>
typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator[](int j){
    return matr_[j];
}

template< typename T>
const typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator[] (int j) const{
    return matr_[j];
}
template<typename T>
typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator()(int i , int j){
    return matr_[i+nrows_*j];
}

template<typename T>
const typename HelperMatrix<T>::value_type& HelperMatrix<T>::operator() (int i, int j) const{
    return matr_[i+nrows_*j];
}

template <typename T>
std::ostream& operator<< (std::ostream& stream_, const HelperMatrix<T>& matr_){
    for(int i = 0; i < matr_.nrows_; ++i){
        for(int j =0; j< matr_.ncols_; ++j ) {
            stream_ << "(" << i << ", " << j << ")" << " :";
            for (const auto elem : matr_(i, j))
                stream_ << elem << " ";
            stream_<<std::endl;
        }
    }
    return stream_;
}

template < typename T>
int HelperMatrix<T>::size(int i )const{
    return matr_[i].size();
}

template < typename T>
int HelperMatrix<T>::size(int i, int j) const{
    return matr_[i + nrows_*j].size();
}

template< typename T>
int HelperMatrix<T>::size() const{
    return nrows_*ncols_;
}

//da rinominare

template <typename T>
std::vector<int> compute_number_elements(const HelperMatrix<T>& matrix_){
    std::vector<int> result(matrix_.size(),0);
    for(int i=0; i<matrix_.size(); ++i)
        result[i] = matrix_.size(i);
    return result;
}

template <typename T>
void helper_neighbors( HelperMatrix<T>& result, const HelperMatrix<T>& node_list){
    for(int node = 0; node < node_list.nrows(); ++node){
        const typename HelperMatrix<T>::value_type& curr_right = node_list(node,0); //questo vector<int> // lista di lati collegati
        // a node da parte "0", i.e da hanno il nodo allora loro dx
        const typename HelperMatrix<T>::value_type& curr_left = node_list(node,1); //questo vector<int> // lista di lati collegati
        // a node da parte "1", i.e da hanno il nodo allora loro sx
        helper_neigh_same_side_(result,curr_right,0);
        helper_neigh_same_side_(result,curr_left,1);

        for(const auto& i : curr_right)
            for(const auto& j : curr_left){
                result(i,0).push_back(j);
                result(j,1).push_back(i);
            }
    }
}

template <typename T>
void helper_neigh_same_side_( HelperMatrix<T>& result, const typename HelperMatrix<T>::value_type& curr,int idx){

    if( !curr.empty() ) {
        for (int i = 0; i < curr.size() - 1; ++i) {
            for (int j = i + 1; j < curr.size(); ++j) {
                result(curr[i], idx).push_back(curr[j]);
                result(curr[j], idx).push_back(curr[i]);
            }
        }
    }
}

#endif //__HELPER_MATRIX_H__