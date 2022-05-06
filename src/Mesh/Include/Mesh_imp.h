#ifndef MESH_IMP_H_
#define MESH_IMP_H_

template <UInt ORDER, UInt mydim, UInt ndim>
MeshHandler<ORDER,mydim,ndim>::MeshHandler(SEXP Rmesh, UInt search) :
	points_(VECTOR_ELT(Rmesh, 0)), sides_(VECTOR_ELT(Rmesh, 6)),
		elements_(VECTOR_ELT(Rmesh, 3)), neighbors_(VECTOR_ELT(Rmesh, 8)),
		 	search_(search) {
		 		if(search==2)
		 			tree_ptr_=make_unique<const ADTree<meshElement> > (Rmesh);
				}

template <UInt ORDER, UInt mydim, UInt ndim>
Point<ndim> MeshHandler<ORDER,mydim,ndim>::getPoint(const UInt id) const
{
	return Point<ndim>(id, points_);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getElement(const UInt id) const
{
	typename meshElement::elementPoints elPoints;
	for (int j=0; j<how_many_nodes(ORDER,mydim); ++j)
		elPoints[j] = getPoint(elements_(id,j));
	return meshElement(id, elPoints);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::getNeighbors(const UInt id_element, const UInt number) const
{
	int id_neighbor{neighbors_(id_element, number)};
	//return empty element if "neighbor" not present (out of boundary!)
	return (id_neighbor==-1) ? meshElement() : getElement(id_neighbor);
}

template <UInt ORDER, UInt mydim, UInt ndim>
template <bool isManifold>
typename std::enable_if<!isManifold, typename MeshHandler<ORDER,mydim,ndim>::meshElement>::type  
MeshHandler<ORDER,mydim,ndim>::findLocation(const Point<ndim>& point) const {
	if(search_==2)
		return findLocationTree(point);
	else if(search_==3)
		return findLocationWalking(point, getElement(0));
	else
		return findLocationNaive(point);
}

template <UInt ORDER, UInt mydim, UInt ndim>
template <bool isManifold>
typename std::enable_if<isManifold, typename MeshHandler<ORDER,mydim,ndim>::meshElement>::type  
MeshHandler<ORDER,mydim,ndim>::findLocation(const Point<ndim>& point) const {
	if(search_==2)
		return findLocationTree(point);
	else
		return findLocationNaive(point);
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::findLocationNaive(const Point<ndim>& point) const
{
	for(UInt id=0; id < elements_.nrows(); ++id){
		meshElement current_element{getElement(id)};
		if(current_element.isPointInside(point))
			return current_element;
	}
	return meshElement(); //default element with NVAL ID
}

// Visibility walk algorithm which uses barycentric coordinate [Sundareswara et al]
//Starting triangles usually n^(1/3) points
template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::findLocationWalking(const Point<ndim>& point, const meshElement& starting_element) const
{
	static_assert(mydim==ndim, "ERROR! WALKING SEARCH CANNOT BE USED ON MANIFOLD MESHES! See mesh_imp.h");
	meshElement current_element{starting_element};
	//Test for found Element, or out of border
	while(current_element.hasValidId() && !current_element.isPointInside(point))
		current_element = getNeighbors(current_element.getId(), current_element.getPointDirection(point));

	return current_element;
}

template <UInt ORDER, UInt mydim, UInt ndim>
typename MeshHandler<ORDER,mydim,ndim>::meshElement MeshHandler<ORDER,mydim,ndim>::findLocationTree(const Point<ndim>& point) const {
	std::set<int> found;
	std::vector<Real> region;
	region.reserve(2*ndim);

	for (UInt i=0; i<ndim; ++i){
		region.push_back(point[i]);
	}
	for (UInt i=0; i<ndim; ++i){
		region.push_back(point[i]);
	}

	if(!tree_ptr_->search(region, found)) {
		return meshElement();
	}

	for (const auto &i : found) {
		const UInt index = tree_ptr_->pointId(i);
		meshElement tmp = getElement(index);
		if(tmp.isPointInside(point)) {
			return tmp;
		}
	}
	return meshElement();
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printPoints(std::ostream& os) const
{
	os<<"# Nodes: "<<points_.nrows()<<std::endl;
	for(UInt i=0; i<points_.nrows(); ++i)
		os<<getPoint(i);
}


template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printElements(std::ostream& os) const
{
	os << "# Triangles: "<< elements_.nrows() <<std::endl;
	for (UInt i = 0; i < elements_.nrows(); ++i )
		os<<getElement(i);
}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printNeighbors(std::ostream& os) const
{
	os << "# Neighbors list: "<< elements_.nrows() <<std::endl;
	for (UInt i = 0; i < elements_.nrows(); ++i ){
		for( UInt j = 0; j < mydim+1; ++j)
			os<<neighbors(i,j)<<" ";
		os<<std::endl;
	}
}

template <UInt ORDER, UInt mydim, UInt ndim>
void MeshHandler<ORDER,mydim,ndim>::printTree(std::ostream & os) const
{
	os << "# Tree characteristic: " <<std::endl;
	if (tree_ptr_)
		os << *tree_ptr_ << std::endl;
	else
		os << "No tree!" <<std::endl;
}

//Linear Networks partial specialization
template<UInt ORDER>
MeshHandler<ORDER,1,2>::MeshHandler(SEXP Rmesh, UInt search) :
        points_(VECTOR_ELT(Rmesh, 0)),
        elements_(VECTOR_ELT(Rmesh, 3)), neighbors_(VECTOR_ELT(Rmesh, 8)),
        search_(search) {

    if(search==2)
        tree_ptr_=make_unique<const ADTree<meshElement> > (Rmesh);

}

template <UInt ORDER>
Point<2> MeshHandler<ORDER,1,2>::getPoint(const UInt id) const
{
    return Point<2>(id, points_);
}

template <UInt ORDER>
typename MeshHandler<ORDER,1,2>::meshElement MeshHandler<ORDER,1,2>::getElement(const UInt id) const
{
    typename meshElement::elementPoints elPoints;
    for (int j=0; j<how_many_nodes(ORDER,1); ++j)
        elPoints[j] = getPoint(elements_(id,j));
    return meshElement(id, elPoints);
}

template <UInt ORDER>
std::vector<typename MeshHandler<ORDER,1,2>::meshElement > MeshHandler<ORDER,1,2>::getNeighbors(const UInt id_element, const UInt number) const
{
    const RIntegerMatrix &id_neighbors = neighbors_(id_element, number);
    std::vector< meshElement > elem_neighbors;
    elem_neighbors.reserve(id_neighbors.nrows()*id_neighbors.ncols());

    for(UInt i=0; i<id_neighbors.nrows()*id_neighbors.ncols(); ++i){
        elem_neighbors.push_back(getElement(id_neighbors[i]));
    }
    //return empty vector if "neighbors" not present (out of boundary!)
    return elem_neighbors;
}

template <UInt ORDER>
typename MeshHandler<ORDER,1,2>::meshElement MeshHandler<ORDER,1,2>::findLocation(const Point<2>& point) const {
    if(search_==2) {
        return findLocationTree(point);
    }
    else
        return findLocationNaive(point);
}

template <UInt ORDER>
typename MeshHandler<ORDER,1,2>::meshElement MeshHandler<ORDER,1,2>::findLocationNaive(const Point<2>& point) const
{
    for(UInt id=0; id < elements_.nrows(); ++id){
        meshElement current_element{getElement(id)};
        if(current_element.isPointInside(point))
            return current_element;
    }
    return meshElement(); //default element with NVAL ID
}

template <UInt ORDER>
typename MeshHandler<ORDER,1,2>::meshElement MeshHandler<ORDER,1,2>::findLocationTree(const Point<2>& point) const {
    std::set<int> found;
    std::vector<Real> region;
    region.reserve(2*2);

    for (UInt i=0; i<2; ++i){
        region.push_back(point[i]);
    }
    for (UInt i=0; i<2; ++i){
        region.push_back(point[i]);
    }

    if(!tree_ptr_->search(region, found)) {
        return meshElement();
    }

    for (const auto &i : found) {
        const UInt index = tree_ptr_->pointId(i);
        meshElement tmp = getElement(index);
        if(tmp.isPointInside(point)) {
            return tmp;
        }
    }
    return meshElement();
}


template <UInt ORDER>
void MeshHandler<ORDER,1,2>::printPoints(std::ostream& os) const
{
    os<<"# Nodes: "<<points_.nrows()<<std::endl;
    for(UInt i=0; i<points_.nrows(); ++i)
        os<<getPoint(i);
}


template <UInt ORDER>
void MeshHandler<ORDER,1,2>::printElements(std::ostream& os) const
{
    os << "# Edges: "<< elements_.nrows() <<std::endl;
    for (UInt i = 0; i < elements_.nrows(); ++i )
        os<<getElement(i);
}

template <UInt ORDER>
void MeshHandler<ORDER,1,2>::printNeighbors(std::ostream& os) const
{
    os << "# Neighbors list: "<< elements_.nrows() <<std::endl;
    for (UInt i = 0; i < elements_.nrows(); ++i){
        for( UInt j = 0; j < 2; ++j){
            for(UInt k=0; k<neighbors_(i,j).nrows()*neighbors_(i,j).ncols(); ++k) {
                os << neighbors_(i, j)[k] << " ";
            }
            os<<"\t";
        }
        os << std::endl;
    }
}

template <UInt ORDER>
void MeshHandler<ORDER,1,2>::printTree(std::ostream & os) const
{
    os << "# Tree characteristic: " <<std::endl;
    if (tree_ptr_)
        os << *tree_ptr_ << std::endl;
    else
        os << "No tree!" <<std::endl;
}


#endif
