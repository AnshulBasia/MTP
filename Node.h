#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string.h>
#include <math.h>

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

using namespace std;

class Node{

public:

	Node(int _NodeId);

	void addEdge(int _destId, bool _direction, int _label);
    void sort_labels();
    void add_label(int l);
	void removeEdge(int _destId, bool _direction, int _label);
    int nodeId;
    vector<int> labels;

	unordered_map < int, int > numFwdEdges;
	unordered_map < int, int > numBwdEdges;

	unordered_map < int, vector < int > > fwd_labelled_edges;
	unordered_map < int, vector < int > > bwd_labelled_edges;
};

Node::Node(int _NodeId){
	nodeId = _NodeId;
}

void Node::addEdge(int _destId, bool _direction, int _label = 0){
	if ( _direction) {
		fwd_labelled_edges[_label].push_back(_destId);
		numFwdEdges[_label] += 1;
	}
	else {
		bwd_labelled_edges[_label].push_back(_destId);
		numBwdEdges[_label] += 1;
	}
}

void Node::removeEdge(int _destId, bool _direction, int _label = 0){
	if ( _direction) {
		fwd_labelled_edges[_label].erase(remove(fwd_labelled_edges[_label].begin(), fwd_labelled_edges[_label].end(), _destId), fwd_labelled_edges[_label].end());
		numFwdEdges[_label] -= 1;
	}
	else {
		bwd_labelled_edges[_label].erase(remove(bwd_labelled_edges[_label].begin(), bwd_labelled_edges[_label].end(), _destId), bwd_labelled_edges[_label].end());
		numBwdEdges[_label] -= 1;
	}
}

void Node::add_label(int l){
	labels.push_back(l);
}

void Node::sort_labels(){
	sort(labels.begin(),labels.end());
}