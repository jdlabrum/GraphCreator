#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

struct vert{
	int id;
	int parent = -1;
	vector<int> edges;
	vector<int> weights;
};

class graph{
public:
	vector<int> nodeIDs;
	vector<vert*> nodes;
	int graphWeight = 0;

	int getNode(int id){
		int i = 0;
		int vsize = nodeIDs.size();
		for (i = 0; i < vsize; i++)
			if (id == nodeIDs[i]) return i;
		return -1;
	}

	int getEdge(vert *vertA, int id){
		int i = 0;
		int verts = vertA->edges.size();
		for (i = 0; i < verts; i++)
			if (id == vertA->edges[i]) return i;
		return -1;
	}

	bool isNode(int id){
		int i = 0;
		int vsize = nodeIDs.size();
		if (vsize == 0) return false;
		else {
			for (i = 0; i < vsize; i++) if (id == nodeIDs[i]) return true;
			return false;
		}
	}

	void addNode(int id){
		if (isNode(id) == false){
				nodeIDs.push_back(id);
				vert *knewVert;
				knewVert = new vert;
				knewVert->id = id;
				nodes.push_back(knewVert);
		}
	}

	void addEdge(int nodeA, int nodeB, int weight){
		int nA = getNode(nodeA);
		int nB = getNode(nodeB);
		if (nA < 0 || nB < 0) cout << "Error adding Edge.";
		else{
			nodes[nA]->edges.push_back(nodeB);
			nodes[nA]->weights.push_back(weight);
			nodes[nB]->edges.push_back(nodeA);
			nodes[nB]->weights.push_back(weight);
			graphWeight += weight;
		}
	}

	void deleteEdge(int nodeA, int nodeB){
		int nA = getNode(nodeA);
		int nB = getNode(nodeB);
		if (nA < 0 || nB < 0) cout << "No such edge.";
		else {
			int edgeA = getEdge(nodes[nA], nodeB);
			int edgeB = getEdge(nodes[nB], nodeA);
			if (edgeA < 0 || edgeB < 0) cout << "No such edge!";
			else{
				nodes[nA]->edges.erase(nodes[nA]->edges.begin() + edgeA);
				nodes[nA]->weights.erase(nodes[nA]->weights.begin() + edgeA);
				nodes[nB]->edges.erase(nodes[nB]->edges.begin() + edgeB);
				nodes[nB]->weights.erase(nodes[nB]->weights.begin() + edgeB);
				cout << "Deleted Edge: " << nodeA << "<->" << nodeB << endl;
			}
		}
	}

	void deleteNode(int id){
		
	//calls deleteEdge for all edges connected to the given id, then deletes node.
		
		if (isNode(id) == true){
			int node = getNode(id);
			int numEdges = nodes[node]->edges.size();
			cout << numEdges << endl;
			for (int i = 0; i < numEdges; i++){
				deleteEdge(id, nodes[node]->edges[0]);
			}
			delete nodes[node];
			nodeIDs.erase(nodeIDs.begin() + node);
			nodes.erase(nodes.begin() + node);
		}
	}

	void printGraph(ofstream& outf){
		
//prints graph with nodes and edges sorted by id
		
		vector<int> tempIDs = nodeIDs;
		vector<int> tempNodes;
		int count;
		sort(tempIDs.begin(), tempIDs.end());
		int nSize = nodeIDs.size();
		for (int i = 0; i < nSize; i++){
			int temp = getNode(tempIDs[i]);
			if (i>0) outf << endl;
			outf << right << setw(2) << tempIDs[i] << ": ";
			tempNodes = nodes[temp]->edges;
			sort(tempNodes.begin(), tempNodes.end());
			count = 0;
			int tnSize = tempNodes.size();
			for (int j = 0; j < tnSize; j++){
				int nID = tempNodes[j];
				for (int k = 0; k < tnSize; k++){
					if (nID == nodes[temp]->edges[k]){
						if (count > 0) outf << ",";
						outf << " " << nID << ":" << nodes[temp]->weights[k];
						count++;
					}
				}
			}
		}
	}
};

graph readGraph(){
	graph newGraph;
	ifstream inf("graphData.txt");
	string input; 
	int edge, nodeA, nodeB;
	while (!inf.eof()){
		inf >> input;
		input = input.substr(1);
		nodeA = atoi(input.c_str());
		newGraph.addNode(nodeA);
		inf >> nodeB;
		newGraph.addNode(nodeB);
		inf >> input;
		input.pop_back();
		edge = atoi(input.c_str());
		newGraph.addEdge(nodeA, nodeB, edge);
	}
		return newGraph;
}

graph minSpanningTree(graph startGraph, int start){

//shortestPath creates a new graph, a shortest path tree, from 
//the passed in' startGraph' rooted at 'start'

	graph mTree;
	vector<int> aNodes;
	vector<int> bNodes;
	vector<int> weights;
	int i = 0;
	int j = 0;
	int t = 0;
	mTree.addNode(start);
	j = startGraph.nodes[startGraph.getNode(start)]->edges.size();

//Starts by adding the 'start' node. 
//all of the start node's edge information is added to the parallel vectors above. 

	for (i = 0; i < j; i++){
		aNodes.push_back(start);
		bNodes.push_back(startGraph.nodes[startGraph.getNode(start)]->edges[i]);
		weights.push_back(startGraph.nodes[startGraph.getNode(start)]->weights[i]);
	}

//Finds the edge with the least weight. Adds the node it connects to and the edge itself.
//Then adds edges from the new node to the parallel vectors as long as they do not connect
//a node that is already in the new graph.

	while (mTree.nodeIDs.size() < startGraph.nodeIDs.size()) {
		j = weights[0];
		t = 0;
		for (i = 0; i < weights.size(); i++){
			if (j > weights[i]){
				t = i;
				j = weights[i];
			}
		}
		if (mTree.isNode(bNodes[t]) == true){
			aNodes.erase(aNodes.begin() + t);
			bNodes.erase(bNodes.begin() + t);
			weights.erase(weights.begin() + t);
		}
		else {
			mTree.addNode(bNodes[t]);
			mTree.addEdge(aNodes[t], bNodes[t], weights[t]);
			for (i = 0; i < startGraph.nodes[startGraph.getNode(bNodes[t])]->edges.size(); i++){
				if (mTree.isNode(startGraph.nodes[startGraph.getNode(bNodes[t])]->edges[i]) == false){
					aNodes.push_back(bNodes[t]);
					bNodes.push_back(startGraph.nodes[startGraph.getNode(bNodes[t])]->edges[i]);
					weights.push_back(startGraph.nodes[startGraph.getNode(bNodes[t])]->weights[i]);
				}
			}
			aNodes.erase(aNodes.begin() + t);
			bNodes.erase(bNodes.begin() + t);
			weights.erase(weights.begin() + t);
		}
	}
	return mTree;
}

int treeSum(graph sGraph){

//Each graph tracks it's total weight. This is updated each time an edge is added. Returns that count.
	
	return sGraph.graphWeight;
}

graph shortestPathTree(graph startGraph, int start){

//shortestPath creates a new graph, a shortest path tree, from 
//the passed in' startGraph' rooted at 'start' 

	graph shortestGraph;
	vector<int> aNodes;
	vector<int> bNodes;
	vector<int> weights;
	vector<int> deletions;
	int i = 0;
	int j = 0;
	int t = 0;
	int check;
	int knew = 0;
	int distance = 0;
	shortestGraph.addNode(start);
	j = startGraph.nodes[startGraph.getNode(start)]->edges.size();

//Starts by adding the 'start' node. 
//all of the start node's edge information is added to the parallel vectors above. 

	for (i = 0; i < j; i++){
		aNodes.push_back(start);
		bNodes.push_back(startGraph.nodes[startGraph.getNode(start)]->edges[i]);
		weights.push_back(startGraph.nodes[startGraph.getNode(start)]->weights[i]);
	}

//finds the edge with the least weight. Adds the node it connects to and the edge itself.

	while (shortestGraph.nodeIDs.size() < startGraph.nodeIDs.size()) {
		j = weights[0];
		t = 0;
		for (i = 0; i < weights.size(); i++){
			if (j > weights[i]){
				t = i;
				j = weights[i];
			}
		}
		if (shortestGraph.isNode(bNodes[t]) == true){
			aNodes.erase(aNodes.begin() + t);
			bNodes.erase(bNodes.begin() + t);
			weights.erase(weights.begin() + t);
		}
		else {
			shortestGraph.addNode(bNodes[t]);
			shortestGraph.addEdge(aNodes[t], bNodes[t], weights[t]);
			shortestGraph.nodes[shortestGraph.getNode(bNodes[t])]->parent = aNodes[t];
			knew = bNodes[t];
			distance = weights[t];
			int j = bNodes.size();

//Next it adds edges from the new node to the parallel vectors.
//If an edge connects to a node not in the the current vectors or graph it adds it as a potential edge. 
//If the edge connects to a node that is in the vectors, it compares the weight of the existing potential
//edge with the weight of the new edge. If the new edge is smaller, the old potential edge is marked for deletion
// and the new potential edge is added. 
		
			for (i = 0; i < startGraph.nodes[startGraph.getNode(knew)]->edges.size(); i++){
				int j = bNodes.size();
				check = -1;
				if (shortestGraph.isNode(startGraph.nodes[startGraph.getNode(knew)]->edges[i]) == false){
					for (int m = 0; m < j; m++){
						if (startGraph.nodes[startGraph.getNode(knew)]->edges[i] == bNodes[m]) check = m;
					}
					if (check < 0){
						aNodes.push_back(knew);
						bNodes.push_back(startGraph.nodes[startGraph.getNode(knew)]->edges[i]);
						weights.push_back(startGraph.nodes[startGraph.getNode(knew)]->weights[i]+distance);
					}
					else {
						if (startGraph.nodes[startGraph.getNode(knew)]->weights[i] + distance < weights[check]){
							aNodes.push_back(knew);
							bNodes.push_back(startGraph.nodes[startGraph.getNode(knew)]->edges[i]);
							weights.push_back(startGraph.nodes[startGraph.getNode(knew)]->weights[i]+distance);
							deletions.push_back(check);
						}
					}
				}

//deletes the edge that was added to the graph from the parallel vectors
//as well as any edges marked for deletion. 

			}
			deletions.push_back(t);
			sort(deletions.rbegin(), deletions.rend());
			for (i = 0; i < deletions.size(); i++){
				aNodes.erase(aNodes.begin() + deletions[i]);
				bNodes.erase(bNodes.begin() + deletions[i]);
				weights.erase(weights.begin() + deletions[i]);
			}
			deletions.clear();
		}
	}
	return shortestGraph;
}

void printPath(graph pathTree, int startNode, int endNode, ofstream& outf){
	vector<int> path;
	int current = endNode;
	int distance = pathTree.nodes[pathTree.getNode(current)]->weights[0];
	
	while (current != startNode){
		path.push_back(current);
		current = pathTree.nodes[pathTree.getNode(current)]->parent;
	}
	
	outf << "Shortest Path from " << startNode << " to " << endNode << ": " << startNode;
	for (int i = 0; i < path.size(); i++) outf << "->" << path[path.size()-(i+1)];
	outf << "\nTotal Path Weight: " << distance << "\n\n";
}

int main(){
	
	ofstream outf("Output.txt");
	graph projectGraph = readGraph();

	outf << "i) Minimum Spanning Tree: \n\n";

	graph minTree = minSpanningTree(projectGraph, 0);
	minTree.printGraph(outf);

	outf << "\n\nTotal Weight: " << treeSum(minTree);
	
	outf << "\n\n\nii) Shortest Path Tree:\n\n";

	graph pathTree = shortestPathTree(projectGraph, 0);
	printPath(pathTree, 0, 26, outf);
	printPath(pathTree, 0, 18, outf);
	printPath(pathTree, 0, 12, outf);

}