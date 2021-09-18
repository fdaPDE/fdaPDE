#ifndef __MESH_INPUT_HELPER_IMP_H__
#define __MESH_INPUT_HELPER_IMP_H__


template<UInt mydim>
template<std::size_t SIZE>
void simplex_container<mydim>::fill_container(const std::array<UInt, SIZE>& ORDERING){
 static_assert(SIZE==mydim*(mydim+1) || (mydim==2 && SIZE==12),
        "ERROR! ORDERING SIZE SHOULD BE EQUAL TO 2X THE NUMBER OF EDGES OR 3X THE NUMBER OF FACES! See: mesh_input_helper_imp.h");

 const UInt num_elements=elements.nrows();

 simplexes.reserve(num_elements*ORDERING.size()/mydim);

 {
   std::array<UInt,mydim> curr;
   for(UInt i=0; i<num_elements; ++i){
     for(UInt j=0; j<ORDERING.size()/mydim; ++j){
        for(UInt k=0; k<mydim; ++k)
          curr[k]=elements(i, ORDERING[mydim*j+k]);
        std::sort(curr.begin(), curr.end());
        simplexes.emplace_back(i,j,curr);
     }
   }
 }

  bin_sort();
  check_duplicates();
  store_indexes();

}

template<UInt mydim>
void simplex_container<mydim>::bin_sort(){

  std::vector<UInt> positions;
  positions.reserve(simplexes.size());
  for(UInt i=0; i<simplexes.size(); ++i)
    positions.push_back(i);

  bin_sort_(mydim-1, positions);

  for(UInt i=0; i<positions.size(); ++i){
    UInt curr=i;
    while(i!=positions[curr]){
      UInt next=positions[curr];
      std::swap(simplexes[curr],simplexes[next]);
      positions[curr]=curr;
      curr=next;
    }
    positions[curr]=curr;
  }
}

// Recursive unction to sort container by ascending #(index+1) element of the arrays
template<UInt mydim>
void simplex_container<mydim>::bin_sort_(const UInt index, std::vector<UInt> &positions){
  // Note the scoping to avoid unnecessary storage!
  {
    std::vector<UInt> offsets{compute_offsets(index, positions)};
    for(UInt i=0; i<positions.size(); ++i){
      while(i!=offsets[i]){
        UInt next=offsets[i];
        std::swap(positions[i],positions[next]);
        std::swap(offsets[i],offsets[next]);
      }
    }
  }

  if(index>0)
    bin_sort_(index-1, positions);
}

template<UInt mydim>
std::vector<UInt> simplex_container<mydim>::compute_offsets(const UInt index, std::vector<UInt> &positions){
  
  const UInt num_points=nodes.nrows();

  std::vector<UInt> counts(num_points, 0);
  for(auto const &pos : positions)
    ++counts[simplexes[pos][index]];

  UInt offset{0};
  for (auto &count : counts){
    UInt curr{count};
    count=offset;
    offset+=curr;
  }


  std::vector<UInt> offsets;
  offsets.reserve(positions.size());
  for (auto const &pos : positions)
    offsets.push_back(counts[simplexes[pos][index]]++);
  return offsets;

}


template<UInt mydim>
void simplex_container<mydim>::check_duplicates(){
  duplicates.reserve(simplexes.size());
  // First face/edge cannot be a duplicate!
  duplicates.push_back(false);
  for(auto it=std::next(simplexes.cbegin()); it!=simplexes.cend(); ++it)
    duplicates.push_back(*std::prev(it) == *it);
}


template<UInt mydim>
void simplex_container<mydim>::store_indexes(){
  distinct_indexes.reserve(std::count(duplicates.begin(), duplicates.end(), false));
  for(UInt i=0; i<duplicates.size(); ++i)
    if(!duplicates[i])
      distinct_indexes.push_back(i);
}



template<UInt mydim>
void simplex_container<mydim>::assemble_subs(SEXP Routput, UInt index) const {
  
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, distinct_indexes.size(), mydim));
  RIntegerMatrix subsimplexes(VECTOR_ELT(Routput, index));

  for(UInt j=0; j<mydim; ++j)
    for(UInt i=0; i<distinct_indexes.size(); ++i)
      subsimplexes(i,j) = simplexes[distinct_indexes[i]][j] + 1;

}

template<UInt mydim>
void simplex_container<mydim>::mark_boundary(SEXP Routput, UInt index) const {
  
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(LGLSXP, distinct_indexes.size(), 1));
  RIntegerMatrix boundarymarkers(VECTOR_ELT(Routput, index));

  for(UInt i=0; i<distinct_indexes.size()-1; ++i)
    boundarymarkers[i] = !duplicates[distinct_indexes[i]+1];

  //Special attention for the last simplex!
  boundarymarkers[distinct_indexes.size()-1] = distinct_indexes.back()+1==duplicates.size() || !duplicates[distinct_indexes.back()+1];
}


template<UInt mydim>
void simplex_container<mydim>::compute_neighbors(SEXP Routput, UInt index) const {
  
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, simplexes.size()/(mydim+1), mydim+1));
  RIntegerMatrix neighbors(VECTOR_ELT(Routput, index));

  for (UInt i=0; i<simplexes.size(); ++i)
    neighbors[i]=-1;

  auto rep_it=duplicates.cbegin();
  simplex_t prev{simplexes.front()};
  for (auto const &curr : simplexes){
    // Note: the first simplex cannot be a duplicate!
    if (*(rep_it++)){
      neighbors(curr.i(), curr.j()) = prev.i()+1;
      neighbors(prev.i(), prev.j()) = curr.i()+1;
    }
    prev=curr;
  }
}

template<UInt mydim>
void simplex_container<mydim>::order2extend(SEXP Routput, UInt index) const {
  static_assert(mydim==2, 
    "ERROR! ORDER 2 EXTENSIONS IS INTENDED FOR EDGE CONTAINERS ONLY! See mesh_input_helper_imp");
  
  const UInt num_extra_nodes = (isTriangleContainer) ? 3 : 6;

  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, simplexes.size()/num_extra_nodes, num_extra_nodes));
  RIntegerMatrix edges_extended(VECTOR_ELT(Routput, index));

  {
    UInt offset{nodes.nrows()};
    UInt pos=0;
    for(auto const &curr : simplexes){
      offset += !duplicates[pos++];
      edges_extended(curr.i(), curr.j()) = offset;
    }
  }
}

//Need specialization, one element could have at most (#edges - 1) neighbors
/*
template<>
void simplex_container<1>::compute_neighbors(SEXP Routput, UInt index) const {
  
  //Upper Triangular part (order 1 not needed ? )
  UInt num_nodes = this->get_num_points();
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, num_nodes * (num_nodes - 1)/2, 1));
  RIntegerMatrix adjency_matrix(VECTOR_ELT(Routput, index));

   
  for (UInt i=0; i<num_nodes*(num_nodes-1)/2; ++i)
    adjency_matrix[i]= 0;

  for(UInt node = 0; node<elements.nrows();++node){
    UInt i = elements(node,0);
    UInt j = elements(node,1);
    UInt k = num_nodes*(num_nodes - 1)/2 - (num_nodes - i) *(num_nodes- i - 1)/2 + j - i - 1;
    adjency_matrix[k] = 1;
    }
  // 
  //computing neighbors matrix (left-neighbors only)
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(INTSXP, elements.size(),elements.size()) );
  RIntegerMatrix neighbors(VECTOR_ELT(Routput, index + 1));
  
  for (UInt i=0; i<simplexes.size(); ++i)
    neighbors[i]=0;
  
  auto rep_it=duplicates.cbegin();
  simplex_t prev{simplexes.front()};
  for (auto const &curr : simplexes){
    // Note: the first simplex cannot be a duplicate!
    if (*(rep_it++)){
      if(curr.j() == 1){ //Forse è un check inutile perché se un nodo è ripetuto
                        //corrisponde sicuramente ad un elemento che sta alla sx
                        //di quello a cui il nodo corrente appartiene 
      neighbors(curr.i(), prev.i()) = 1; //elemento curr.i() ha alla sua sx l'elemento prev.i()
      }
    }
    prev=curr;
  }
}
*/


template<>
void simplex_container<1>::compute_neighbors(SEXP Routput, UInt index) const {
  
  //SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(VECSXP, elements.nrows(), 2));

  //     * - - - - *        edge
  //     0         1        nodes
  //    (1)       (0)       sides
  //each row contains two lists of the edges to which the node belongs
  //j = 0 there is the "left" list (current node is on the left of the edge)
  //j = 1 there is the "right" list(current node is on the right of the edge)
  HelperMatrix<int> nodes_neighbors(nodes.nrows(),2);

  HelperMatrix<int> neighbors(elements.nrows(),2);

  for(const auto& curr : simplexes){
      nodes_neighbors(curr[0],curr.j()).push_back( curr.i() ); //Brutto
  }
  //filling neighbors matrix
  helper_neighbors( neighbors, nodes_neighbors);

  auto tmp = compute_number_elements(neighbors);
  RIntegerMatrix number_neighbors(tmp.data(), elements.nrows(), 2);

  //filling R Data Structure
  SET_VECTOR_ELT(Routput, index, Rf_allocMatrix(VECSXP, elements.nrows(), 2));
  for(int i=0; i<2*elements.nrows(); ++i){
      SET_VECTOR_ELT( VECTOR_ELT(Routput, index), i , Rf_allocMatrix(INTSXP, 1, number_neighbors[i]));
      RIntegerMatrix neigh( VECTOR_ELT( VECTOR_ELT(Routput, index), i ));
      for(int j=0; j<number_neighbors[i]; ++j) //BRUTTO, NB "+1" perché numerazione in R parte da 1
          neigh[j] = neighbors[i][j] + 1;

    }
}



#endif
