#include "Node.h"

// get_clock method for time keeping
void get_time_clock(struct timespec * ts) {
	#ifdef __MACH__
		clock_serv_t cclock; 
		mach_timespec_t mts; 
		host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock); 
		clock_get_time(cclock, &mts); 
		mach_port_deallocate(mach_task_self(), cclock); 
		ts->tv_sec = mts.tv_sec; 
		ts->tv_nsec = mts.tv_nsec; 
	#else
	clock_gettime( CLOCK_MONOTONIC, ts); 
	#endif
}

// Main Class Grpaph
class Graph {

public:
	float c_walkLength, c_numWalks; 

	Graph(char *graphfilename, int dir_control); 

	void readGraph(char *graphfilename); 
	void readFeat(char *featfilename); 

	void readAndRunQuery(char *queryFileName); 

	void addEdge(int src, int dst, int label, int dir_control); 
	int diameter(); 

	int rr(int src, int dst, int query_type, vector<int> query_labels); 
	int bbfs(int src, int dst, int query_type, vector<int> query_labels); 
	double charactersticPathLength();
	int floydWarshall();
	void Greedy(int k);
	void ModifiedGreedy(int k);
	void removeEdge(int m, int n);
	
private:

	struct timespec start, finish;
	double initializationTime; 

	int rand_index; 
	vector<int> rand_vector; 

	vector<Node*> nodes; 
	int numNodes; 
	int numEdges; 
	int **dist;
	int numQueries;

	bool rrParamsInitialized; 
	int numWalks, numStops, walkLength; 

}; 

Graph::Graph(char *graphfilename, int dir_control) {

	string line; 
	initializationTime = 0; 

	double elapsed;

	get_time_clock( &start); 

    numNodes = 0; 

	
	// Add edges between Nodes, directed or undirected controlled 
	ifstream myfile2(graphfilename); 
	while (getline(myfile2, line)) {

		if (line[0] == '#') continue; 
		char *cstr = &line[0u]; 

		char *t = strtok(cstr, " "); 
		int u = atoi(t); 
		t = strtok(NULL, " "); 
		int v = atoi(t); 
		//if(u>=10 || v>=10)continue;
		if(u > 200 && u < 300 && v > 200 && v < 300){
			u %= 100;
			v %= 100;
		}
		else continue;
		for (int j = numNodes; j <= u; j ++ ) {
            nodes.push_back(new Node(j)); 
			numNodes ++ ; 
        }

		

		for (int j = numNodes; j <= v; j ++ ) {
            nodes.push_back(new Node(j)); 
			numNodes ++ ; 
        }

		int l; 
		t = strtok(NULL, " "); 
		if (t == NULL) {
			l = 0; 
		}
		else 
			l = atoi(t); 

		addEdge(u, v, l, dir_control); 
	}
	myfile2.close();

	// Create random number array
	rand_vector.reserve(numNodes);
	for (int i = 0; i < numNodes; i ++) {
		rand_vector.push_back(rand() % numNodes); 
	}
	rand_index = 0; 

	// Estimate diameter of graph by sampling 10 random nodes
	int dia = diameter();
	for (int i = 0; i < 4; i ++ ) {
		int d = diameter(); 
		if (d > dia) dia = d; 
	}

	cout << "NumWalks: " << numWalks << ", WalkLength: " << walkLength << endl; 

	get_time_clock(&finish); 
	initializationTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10, 9); 

	numQueries = 0;
	cout << "1. Initialization Time = " << initializationTime << endl; 

}

// Simple bfs to estimate the depth
int Graph::diameter() {
	vector<int> color;
	color.reserve(numNodes);
	for (int i = 0; i < numNodes; i++) color.push_back(0);
	int src = rand()%numNodes;
	color[src] = 1;
	int dia = 1;
	vector<int> queue;
	int u = src;
	while (u >= 0){
		Node *n = nodes[u];
		int numE = n->numFwdEdges[0];
        vector<int> fwdedges = n -> fwd_labelled_edges[0];
		for (int i = 0; i < numE; i++){
			int v =fwdedges[i];
			if (color[v] >= 1) continue;
			color[v] = color[u] + 1;
			queue.push_back(v);
		}
		if (queue.size() == 0) break;
		u = queue[0];
		queue.erase(queue.begin());
	}
	for (int i = 0; i < numNodes; i++)
		if (dia < (color[i]))
			dia = color[i];
	return dia;
}

// Adding edges constrained on whether directional edge or not
void Graph::addEdge(int src, int dst, int l, int dir_control = 0) {
	if ((src >= numNodes) || (dst >= numNodes)) return; 
	nodes[src]->addEdge(dst, true, l); 
	nodes[dst]->addEdge(src, false, l); 
	if (dir_control == 0) {
		nodes[src]->addEdge(dst, false, l); 
		nodes[dst]->addEdge(src, true, l); 
	}
	numEdges ++ ; 
}

void Graph::readAndRunQuery(char *queryFileName) {

	int total_queries = 0; 
	int tp = 0; 
	int fp = 0; 
	int tn = 0; 
	int fn = 0; 
	int ques = 0; 
 
	double mean_BFS_time = 0, mean_ques_BFS_time = 0, mean_True_BFS_time = 0, mean_False_BFS_time = 0; 
	double mean_speedup = 0, mean_True_speedup = 0, mean_False_speedup = 0; 
	vector<double> BFS_times; 

	double elapsed_bbfs;

	ifstream myfile(queryFileName); 
	string line; 

	// query parsing
	while (getline(myfile, line)) {

		total_queries ++; 

		char *cstr = &line[0u];
		char *t = strtok(cstr, " "); 
		int u = atoi(t); 

		t = strtok(NULL, " "); 
		int v = atoi(t); 

		t = strtok(NULL, " "); 
		int query_type = atoi(t); 

		t = strtok(NULL, " "); 
		int num_labels = atoi(t); 

		vector<int> query_labels; 
		while(num_labels--) {
			t = strtok(NULL, " "); 
			query_labels.push_back(atoi(t)); 
		}

		// bbfs calls
		int x1 = 0; 
        get_time_clock( &start); 
        x1 = bbfs(u, v, query_type, query_labels); 
        get_time_clock( &finish); 

        elapsed_bbfs = (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10, 9); 
		cout << "BBFS: " << x1 << " BBFSTime: " << elapsed_bbfs; 
	}
}

int Graph::bbfs(int src, int dst, int query_type, vector<int> query_labels) {

	// Contains whatever Node has been visited to reach in the current walk
	unordered_set < int > soFar; 

	// Queue contains: Node, label, set of soFar nodes
	vector < pair < int, pair < int, unordered_set<int> > > > queue_F; 
	vector < pair < int, pair < int, unordered_set<int> > > > queue_B; 

	// Current Node -> label number -> vector of paths(vector) to node
	unordered_map < int, unordered_map < int, vector < unordered_set < int > > > > fwd; 
	unordered_map < int, unordered_map < int, vector < unordered_set < int > > > > bwd; 

	// Initialize the queues
	queue_F.push_back(make_pair(src, make_pair(-1, soFar))); 
	queue_B.push_back(make_pair(dst, make_pair(query_labels.size(), soFar))); 

	// If any of the queues becomes empty we have not reachable condition
	while(queue_F.size() > 0 && queue_B.size()>0) {

		// Time limit exceeded
		get_time_clock( &finish); 
		if ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10, 9) > 60) return 2; 

		// Forward direction basic initialization: u->Node, l-> label, soFar-> Nodes already reached
		int u = queue_F[0].first; 
		int l = queue_F[0].second.first; 
		soFar = queue_F[0].second.second; 

		// Pop the first element
		queue_F.erase(queue_F.begin()); 

		// Get where the edges from the Node in the given direction
		vector<int> edges = nodes[u]->fwd_labelled_edges[0]; 
		vector<int> labels_no; 

		// Insert into the Nodes visited so far the current Node
		soFar.insert(u); 

		// change label according to query type
		if (query_type == 2) {
			for (int i = 0; i < query_labels.size(); i ++ ) {
				labels_no.push_back(i); 
			}
		}

		else if (query_type == 3) {
			l ++ ; 
			if (l == query_labels.size())
				l = 0; 
			labels_no.push_back(l); 
		}

		else if (query_type == 4) {
			if (l == -1) {
				labels_no.push_back(0); 
			}
			else {
				labels_no.push_back(l); 
				if (l != query_labels.size() - 1) {
					labels_no.push_back(l + 1); 
				}
			}
		}

		// look at all the possible edges out of the node
		for (int i = 0; i<edges.size(); i ++ ) {
			int e = edges[i]; 

			// Condition path has a cycle
			if (soFar.find(e) != soFar.end()) {
				continue; 
			}

			// Search for query label in the labels of node under question
			vector<int> labels = nodes[e]->labels; 
			for (int i = 0; i < labels_no.size(); i ++ ) {
				l = labels_no[i]; 
				if (binary_search(labels.begin(), labels.end(), query_labels[l])) {

					// No meeting with backward BFS
					if (bwd.find(e) == bwd.end()) {
						fwd[e][l].push_back(soFar); 
						queue_F.push_back(make_pair(e, make_pair(l, soFar))); 
					}

					// Meeting with backward BFS
					else {

						// If same label then do a part by part matching to check whether 
						if (bwd[e].find(l) != bwd[e].end()) {

							// Set to check in
							vector < unordered_set < int > > matches = bwd[e][l]; 
							for (int i = 0; i < matches.size(); i ++) {

								//Iterate over the current set and find in the set to check in
								unordered_set < int > ::iterator it; 
								bool flag = true; 
								for (it = soFar.begin(); it != soFar.end(); it ++ ) {
									if (matches[i].find(*it) != matches[i].end()) {
										flag = false; 
										break; 
									}
								}
								if (flag) {
									return 1; 
								}
							}
						}
						
						fwd[e][l].push_back(soFar); 
						queue_F.push_back(make_pair(e, make_pair(l, soFar))); 
					}
				}
			}
		}

		labels_no.clear(); 

		// Backward walk initialization
		int v = queue_B[0].first; 
		l = queue_B[0].second.first; 
		soFar = queue_B[0].second.second; 

		// Pop the first element
		queue_B.erase(queue_B.begin()); 

		// Get where the edges from the Node in the given direction
		edges = nodes[v]->bwd_labelled_edges[0]; 

		// Insert into the Nodes visited so far the current Node
		soFar.insert(v); 

		// change label according to query type
		if (query_type == 2) {
			for (int i = 0; i < query_labels.size(); i ++ ) {
				labels_no.push_back(i); 
			}
		}

		else if (query_type == 3) {
			l--; 
			if (l == -1)
				l = query_labels.size() - 1; 
			labels_no.push_back(l); 
		}
		
		else if (query_type == 4) {
			if (l == query_labels.size()) {
				labels_no.push_back(l-1); 
			}
			else {
				labels_no.push_back(l); 
				if (l != 0) {
					labels_no.push_back(l - 1); 
				}
			}
		}

		// look at all the possible edges out of the node
		for (int i = 0; i < edges.size(); i ++ ) {
			int e = edges[i]; 

			// Condition path has a cycle
			if (soFar.find(e) != soFar.end()) {
				continue; 
			}

			// Search for query label in the labels of node under question
			vector<int> labels = nodes[e]->labels; 
			for (int i = 0; i < labels_no.size(); i ++ ) {
				l = labels_no[i]; 
				if (binary_search(labels.begin(), labels.end(), query_labels[l])) {

					// No meeting with forward BFS
					if (fwd.find(e) == fwd.end()) {
						bwd[e][l].push_back(soFar); 
						queue_B.push_back(make_pair(e, make_pair(l, soFar))); 
					}

					// Meeting with forward BFS
					else {
						// If same label then do a part by part matching to check whether 
						if (fwd[e].find(l) != fwd[e].end()) {

							// Set to check in
							vector < unordered_set < int > > matches = fwd[e][l]; 
							for (int i = 0; i < matches.size(); i ++ ) {

								//Iterate over the current set and find in the set to check in
								unordered_set < int > ::iterator it; 
								bool flag = true; 
								for (it = soFar.begin(); it != soFar.end(); it ++ ) {
									if (matches[i].find(*it) != matches[i].end()) {
										flag = false; 
										break; 
									}
								}
								if (flag) {
									return 1; 
								}
							}
						}
						
						bwd[e][l].push_back(soFar); 
						queue_B.push_back(make_pair(e, make_pair(l, soFar))); 
					}
				}
			}
		}
	}	
	return 0; 
}

double Graph::charactersticPathLength(){
	cout<<"Computing Characterstic Path Length"<<endl;
	cout<<"Number of Nodes: "<<numNodes<<endl;
	int total = 0;
	vector<int> nu;
	for(int i = 0; i < numNodes; i++){
		if(i%500==0)cout<<i<<endl;
		for(int j = i+1; j < numNodes; j++){
			total = total + bbfs(i,j,2,nu);
		}
	}

	return total;//*2/(numNodes*(numNodes-1));
}

int Graph::floydWarshall() { 
	int V = numNodes;
    /* dist[][] will be the output matrix that will finally have the shortest  
      distances between every pair of vertices */
    // int dist[V][V]
	int i, j, k; 
	// memset(dist, 10000, sizeof(dist[0][0]*V*V));
	
	dist = new int*[V];
	for(int i=0;i<V;i++){
		dist[i] = new int[V];
		for(int j=0;j<V;j++){
			dist[i][j] = 10000;
		}
	}
    /* Initialize the solution matrix same as input graph matrix. Or  
       we can say the initial values of shortest distances are based 
       on shortest paths considering no intermediate vertex. */
	for(int i=0;i<nodes.size();i++){
		int id = nodes[i]->nodeId;
		for(int j=0;j<nodes[i]->fwd_labelled_edges[0].size();j++){
			dist[id][nodes[i]->fwd_labelled_edges[0][j]] = 1;
		}
	}
    // for (i = 0; i < V; i++) 
    //     for (j = 0; j < V; j++) 
    //         dist[i][j] = graph[i][j]; 
  
    /* Add all vertices one by one to the set of intermediate vertices. 
      ---> Before start of an iteration, we have shortest distances between all 
      pairs of vertices such that the shortest distances consider only the 
      vertices in set {0, 1, 2, .. k-1} as intermediate vertices. 
      ----> After the end of an iteration, vertex no. k is added to the set of 
      intermediate vertices and the set becomes {0, 1, 2, .. k} */
    for (k = 0; k < V; k++) 
    { 
        // Pick all vertices as source one by one 
        for (i = 0; i < V; i++) 
        { 
            // Pick all vertices as destination for the 
            // above picked source 
            for (j = 0; j < V; j++) 
            { 
                // If vertex k is on the shortest path from 
                // i to j, then update the value of dist[i][j] 
                if (dist[i][k] + dist[k][j] < dist[i][j]){ 
                    dist[i][j] = dist[i][k] + dist[k][j];
					dist[j][i] = dist[i][j];
				} 
            } 
        } 
    } 
	int total = 0;
    for (int i = 0; i < V; i++) 
    { 
        for (int j = 0; j < V; j++) 
        { 
            //if(dist[i][j]>1)cout<<i<<" "<<j<<endl;
			if (dist[i][j] < 10000) 
                total += dist[i][j]; 
        } 
        
    } 
	return total;
} 

void Graph::Greedy(int k){
	int first,second;
	cout<<numEdges<<endl;
	for(int i=0;i<k;i++){
		cout<<"K = "<<i+1<<endl;
		int temp = floydWarshall();
		cout<<"edges: "<<numEdges<<endl;
		cout<<temp<<endl;
		int orig = temp;
		for(int m=0;m<numNodes;m++){
			//if(m%500==0)cout<<"done with node "<<m<<endl;
			for(int n=m+1;n<numNodes;n++){
				//cout<<"done with node-- "<<n<<endl;
				if(dist[m][n] > 1){
					
					
					// for(int q=0;q<numNodes;q++){
					// 	for(int w=0;w<numNodes;w++){
							
					// 		cout<<dist[q][w]<<" ";
					// 	}
					// 	cout<<endl;
					// }
					
					
					addEdge(m,n,0);
					//cout<<"edges: "<<numEdges<<endl;
					int cl = floydWarshall();
					
					if(cl<temp){
						temp = cl;
						first = m;
						second = n;
					}
					removeEdge(m,n);
				}
			}
		}
		if(orig>temp){
			cout<<first<<" "<<second<<endl;
			addEdge(first,second,0);
		}
		
	}
	cout<<"Final"<<endl;
	int temp = floydWarshall();
	cout<<"edges: "<<numEdges<<endl;
	cout<<temp<<endl;
}

void Graph::ModifiedGreedy(int k){
	int first,second,first2,second2;
	cout<<numEdges<<endl;
	for(int i=0;i<k/2;i++){
		cout<<"K = "<<i+1<<endl;
		int temp = floydWarshall();
		int temp2 = temp;
		cout<<"edges: "<<numEdges<<endl;
		cout<<temp<<endl;
		int orig = temp;
		for(int m=0;m<numNodes;m++){
			//if(m%500==0)cout<<"done with node "<<m<<endl;
			for(int n=m+1;n<numNodes;n++){
				//cout<<"done with node-- "<<n<<endl;
				if(dist[m][n] > 1){
					
					addEdge(m,n,0);
					
					int cl = floydWarshall();
					
					if(cl<temp){
						temp2 = temp;
						first2 = first;
						second2 = second;
						temp = cl;
						first = m;
						second = n;
					}
					else if(cl < temp2){
						temp2 = cl;
						first2 = m;
						second2 = n;
					}
					removeEdge(m,n);
				}
			}
		}
		if(orig>temp){
			cout<<"first min ="<<temp<<endl;
			cout<<first<<" "<<second<<endl;
			addEdge(first,second,0);
		}
		if(orig>temp2){
			cout<<"second min ="<<temp2<<endl;
			cout<<first2<<" "<<second2<<endl;
			addEdge(first2,second2,0);
		}
		
	}
	cout<<"Final"<<endl;
	int temp = floydWarshall();
	cout<<"edges: "<<numEdges<<endl;
	cout<<temp<<endl;
}

void Graph::removeEdge(int src, int dst){
	if ((src >= numNodes) || (dst >= numNodes)) return; 
	nodes[src]->removeEdge(dst, true); 
	nodes[dst]->removeEdge(src, false); 
	
	nodes[src]->removeEdge(dst, false); 
	nodes[dst]->removeEdge(src, true); 
	
	numEdges -= 1 ; 
}

