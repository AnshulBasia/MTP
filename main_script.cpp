#include "Graph.h"

int main(int argc, char *argv[]){

	
	char *graphFile = argv[1];
	
	
	
	srand(time(NULL));
	Graph *tg = new Graph(graphFile, 0);
	//cout<<tg->floydWarshall()<<endl;
	tg->ModifiedGreedy(10);
	
}
