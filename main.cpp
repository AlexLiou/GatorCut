#include <iostream>
#include <ctime>
#include <algorithm>
#include <limits>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include<iostream>
#include<cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include<vector>
#include <sstream>
#include <cmath>
#include <stdlib.h> //needed for atof -> convert string to double
#include "main.h"
#include <iterator>


using namespace std;

LinkedNode::LinkedNode(double distance, Node* neighbor, LinkedNode* next){
    this->distance = distance;
    this->neighbor = neighbor;
    this->next = next;
}

double LinkedNode::getDistance(){
    return this->distance;
}

Node* LinkedNode::getNeighbor(){
    return this->neighbor;
}

LinkedNode* LinkedNode::getNext(){
    return this->next;
}

void LinkedNode::setNext(LinkedNode* next){
    this->next = next;
}

Node::Node(double lon, double lat, string name, vector <Node*> neighbors, vector <double> distances){
    this->lon=lon;
    this->lat=lat;
    this->name = name;
    this->neighbors=neighbors;
    this->distances=distances;
}

string Node::getName(){
    return this->name;
}

double Node::getLon(){
    return this->lon;
}

double Node::getLat(){
    return this->lat;
}

vector <Node*> Node::getNeighbors(){
    return this->neighbors;
}

void Node::pushNeighbors(Node* neighbor){
    this->neighbors.push_back(neighbor);
}

vector <double> Node::getDistances(){
    return this->distances;
}

void Node::pushDistances(double distance){
    this->distances.push_back(distance);
}

vector<string> createTokens(string s){
//    std::cout<<s<<endl;
    vector <string> lineVector;
    string temp="";
    for (int i=0; i<s.length();i++){
        char c = s[i];
        //       std::cout<<c<<" "<<endl;
//        if(c==' '){
//            continue;
//        }
        if(c != ','){
            temp+=c;
            //std::cout<<temp<<endl;
        }
        else{
//            if(temp != " " && !temp.empty()){
//                std::cout<<"THISIS"<<temp<<endl;
            lineVector.push_back(temp);
//            }
            //overload operator
            //I commented out because I think its easier/more efficent without the commas
            /*string jk;
            stringstream ss;
            ss << c;
            ss >> jk;
            lineVector.push_back(jk);*/
            temp = "";
        }
    }
    return lineVector;
}

//checks to see if a string matches the name of a node in the locations vector, changed the value of index to the index where a match is found
bool locationNameExists(string checkCase, vector <Node*> locations, int & index){
    for(int i=0; i<locations.size(); i++){
        if(checkCase==locations[i]->getName())
            index = i;
        return true;
    }
    return false;
}
//checks to see if a string matches the name of a node in the locations vector
bool locationNameExists(string checkCase, vector <Node*> locations){
    for(int i=0; i<locations.size(); i++){
        if(checkCase==locations[i]->getName())
            return true;
    }
    return false;
}

// A structure to represent a node in adjacency list
struct AdjListNode
{
    int dest;
    int weight;
    struct AdjListNode* next;
};

// A structure to represent an adjacency liat
struct AdjList
{
    struct AdjListNode *head;  // pointer to head node of list
};

// A structure to represent a graph. A graph is an array of adjacency lists.
// Size of array will be V (number of vertices in graph)
struct Graph
{
    int V;
    struct AdjList* array;
};

// A utility function to create a new adjacency list node
struct AdjListNode* newAdjListNode(int dest, int weight)
{
    struct AdjListNode* newNode =
            (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}

// A utility function that creates a graph of V vertices
struct Graph* createGraph(int V)
{
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;

    // Create an array of adjacency lists.  Size of array will be V
    graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));

    // Initialize each adjacency list as empty by making head as NULL
    for (int i = 0; i < V; ++i)
        graph->array[i].head = NULL;

    return graph;
}

// Adds an edge to an undirected graph
void addEdge(struct Graph* graph, int src, int dest, int weight)
{
    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    struct AdjListNode* newNode = newAdjListNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;

    // Since graph is undirected, add an edge from dest to src also
    newNode = newAdjListNode(src, weight);
    newNode->next = graph->array[dest].head;
    graph->array[dest].head = newNode;
}

// Structure to represent a min heap node
struct MinHeapNode
{
    int  v;
    int dist;
};

// Structure to represent a min heap
struct MinHeap
{
    int size;      // Number of heap nodes present currently
    int capacity;  // Capacity of min heap
    int *pos;     // This is needed for decreaseKey()
    struct MinHeapNode **array;
};

// A utility function to create a new Min Heap Node
struct MinHeapNode* newMinHeapNode(int v, int dist)
{
    struct MinHeapNode* minHeapNode =
            (struct MinHeapNode*) malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}

// A utility function to create a Min Heap
struct MinHeap* createMinHeap(int capacity)
{
    struct MinHeap* minHeap =
            (struct MinHeap*) malloc(sizeof(struct MinHeap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array =
            (struct MinHeapNode**) malloc(capacity * sizeof(struct MinHeapNode*));
    return minHeap;
}

// A utility function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b)
{
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct MinHeap* minHeap, int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size &&
        minHeap->array[left]->dist < minHeap->array[smallest]->dist )
        smallest = left;

    if (right < minHeap->size &&
        minHeap->array[right]->dist < minHeap->array[smallest]->dist )
        smallest = right;

    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        MinHeapNode *smallestNode = minHeap->array[smallest];
        MinHeapNode *idxNode = minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

// A utility function to check if the given minHeap is empty or not
int isEmpty(struct MinHeap* minHeap)
{
    return minHeap->size == 0;
}

// Standard function to extract minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap* minHeap)
{
    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode* root = minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Function to decreasy dist value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap* minHeap, int v, int dist)
{
    // Get the index of v in  heap array
    int i = minHeap->pos[v];

    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;

    // Travel up while the complete tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist)
    {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i-1)/2;
        minHeap->pos[minHeap->array[(i-1)/2]->v] = i;
        swapMinHeapNode(&minHeap->array[i],  &minHeap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct MinHeap *minHeap, int v)
{
    if (minHeap->pos[v] < minHeap->size)
        return true;
    return false;
}

void printPath(int parent[], int src, int destination, map<int, string> hashmap){
    if(destination==src) return;
    printPath(parent, src, parent[destination], hashmap);
    cout<<hashmap.at(destination)<<"\t\t";
    destination=parent[destination];
}

// A utility function used to print the solution
void printArr(int dist[], int src, int parent[], int destination, map<int, string> hashmap)
{
    printf("Source  Destination   Distance from Source   Path\n");
    cout<<hashmap.at(src)<<"\t\t"<<hashmap.at(destination)<<"\t\t"<<dist[destination]<<"\t\t"<<hashmap.at(src)<<"\t\t";
    printPath(parent, src, destination, hashmap);

}

// The main function that calulates distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
void dijkstra(struct Graph* graph, int src,int dest, map<int, string> hashmap)
{   //starting
    int V = graph->V;// Get the number of vertices in graph
    int dist[V];      // dist values used to pick minimum weight edge in cut
    int destinIndex= dest;
    int parent [V];
    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = INT_MAX;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // Make dist value of src vertex as 0 so that it is extracted first
    minHeap->array[src] = newMinHeapNode(src, dist[src]);
    minHeap->pos[src]   = src;
    dist[src] = 0;
    parent[src]=-1;
    decreaseKey(minHeap, src, dist[src]);

    // Initially size of min heap is equal to V
    minHeap->size = V;

    // In the following loop, min heap contains all nodes
    // whose shortest distance is not yet finalized.
    while (!isEmpty(minHeap))
    {
        // Extract the vertex with minimum distance value
        struct MinHeapNode* minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v; // Store the extracted vertex number

        // Traverse through all adjacent vertices of u (the extracted
        // vertex) and update their distance values
        struct AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl != NULL)
        {
            int v = pCrawl->dest;

            // If shortest distance to v is not finalized yet, and distance to v
            // through u is less than its previously calculated distance
            if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX &&
                pCrawl->weight + dist[u] < dist[v])
            {
                parent[v]=u;
                dist[v] = dist[u] + pCrawl->weight;

                // update distance value in min heap also
                decreaseKey(minHeap, v, dist[v]);
            }
            pCrawl = pCrawl->next;
        }
    }

    // print the calculated shortest distances
    printArr(dist, src, parent, dest, hashmap);
}



int main() {

    // create the graph given in above fugure
    int V = 9;

    struct Graph* graph = createGraph(V);
    addEdge(graph, 0, 1, 4);
    addEdge(graph, 0, 7, 8);
    addEdge(graph, 1, 2, 8);
    addEdge(graph, 1, 7, 11);
    addEdge(graph, 2, 3, 7);
    addEdge(graph, 2, 8, 2);
    addEdge(graph, 2, 5, 4);
    addEdge(graph, 3, 4, 9);
    addEdge(graph, 3, 5, 14);
    addEdge(graph, 4, 5, 10);
    addEdge(graph, 5, 6, 2);
    addEdge(graph, 6, 7, 1);
    addEdge(graph, 6, 8, 6);
    addEdge(graph, 7, 8, 7);

    //dijkstra(graph, 0,2);


    vector <Node*> locations;//vector of node locations
    vector <string> temp;
    vector <string> token;
    vector <vector <int> > adjMatrix;
    vector <LinkedNode*> adjList;
    vector <Node*> neighborsPassToNode;
    vector <double> distancesPassToNode;



map<int,string> hashmap;
    hashmap[0]="StudentRec";
    hashmap[1]="Marston";
    hashmap[2]="Turlington";
    hashmap[3]="LibWest";
    hashmap[4]="CenturyTower";
    hashmap[5]="Potato";
    hashmap[6]="SouthwestGym";
    hashmap[7]="NEB";
    hashmap[8]="CISE";
    std::cout<<"\n"<<hashmap.at(3)<<endl;

map<string, int> maphash;
        maphash["StudentRec"]=0;
        maphash["Marston"]=1;
        maphash["Turlington"]=2;
        maphash["LibWest"]=3;
        maphash["CenturyTower"]=4;
        maphash["Potato"]=5;
        maphash["SouthwestGym"]=6;
        maphash["NEB"]=7;
        maphash["CISE"]=8;



    /*int lines = 0;
    string filename = "map.txt";

    ifstream infile(filename);
    //    infile.open(filename);
    string str;
    while(getline(infile,str)){
        temp = createTokens(str);
        lines++;
        for(int i=0; i<temp.size();i++){
            //if (temp[i] != ","){
            token.push_back(temp[i]);
            // }
        }
    }

    //tests token vector
    /*for(auto p: token){
        std::cout<<p<<endl;
    }
    //iterates through the token list and creates a Node for each one specified
    string name;
    double lon;
    double lat;
    int iter = 0;
    //each line correlates to a new node so each loop creates a node
    for(int i = 1; i <= lines; i++){
        name = token[iter];
        iter++;
        lon = atof(token[iter].c_str()); //converts string to double
        iter++;
        lat= atof(token[iter].c_str());
        if(i!=lines)
            iter++;
        locations.push_back(new Node(lat, lon, name, neighborsPassToNode, distancesPassToNode));
    }
    //prints out to test the Node locations vector
    /*for(int i= 0; i<locations.size(); i++){
    	cout<< locations[i]->getlon();
    	cout<< locations[i]->getY();
    	cout<< locations[i]->getName();
    }

    string adjacencyFileName = "adj.txt";
    //clears token vector for reuse
    token.clear();
    ifstream afile(adjacencyFileName);
    while(getline(afile,str)){
        temp = createTokens(str);
        for(int i=0; i<temp.size();i++){
            //if (temp[i] != ","){
            token.push_back(temp[i]);
            // }
        }
    }


    vector <int> nodeAdj;
    iter = 0;
    //iterates through the token list, knowing each line should have a number of 1s/0s/-1s equal to the locations.size()
    //which is equal to the number of lines in map.txt
    //pushes the 1D vector that corresponds to a node into the 2D vector that holds adjcacency for all nodes
    for(int i = 1; i <= lines; i++){
        for(int j = 0; j<lines; j++){
            nodeAdj.push_back(atoi(token[iter].c_str()));
            iter++;
        }
        adjMatrix.push_back(nodeAdj);
        nodeAdj.clear();
    }
    //tests adjMatrix
    /*for(int x=0;x<adjMatrix.size();x++)
    {
        for(int y=0;y<adjMatrix[x].size();y++)
        {
            cout<<adjMatrix[x][y];
        }
    cout<<endl;
    }

    for(int i = 0; i < lines; i++){

        for(int j = 0; j<lines; j++){
            if(adjMatrix[i][j]==1){
                locations[i]->pushNeighbors(locations[j]);
                locations[i]->pushDistances(getDist(locations[i]->getLat(), locations[i]->getLon(), locations[j]->getLat(), locations[j]->getLon()));
            }

        }

    }

    LinkedNode* emptyLink;
    LinkedNode* tempLink;

    for(int i = 0; i < locations.size(); i++){
        for (int j = 0; j < locations[i]->getNeighbors().size(); ++j)
        {
            if(j==0){
                adjList.push_back(new LinkedNode(locations[i]->getDistances()[j], locations[i]->getNeighbors()[j], emptyLink));
                tempLink = adjList[i];
            }
            else{
                tempLink->setNext(new LinkedNode(locations[i]->getDistances()[j], locations[i]->getNeighbors()[j], emptyLink));
            }
        }
    }*/




    int choice;
    while(choice!=6){
        cout<<"\n--------------------Choices--------------------\n";
        cout<<"                 1. Directions\n";
        cout<<"                2. Add Location\n";
        cout<<"               3. Print Locations\n";
        cout<<"                5. Database test\n";
        cout<<"                    6. Exit\n\n";
        cin >> choice;
        bool fail = cin.fail();
        if(choice > 6 || choice < 1)
            fail = true;
        //validates the user input for being both an integer and in the range of valid input
        while(fail){
            cout<<"Invalid input, please input a whole number between 1 and 6.\n";
            cin.clear();
            cin.ignore(10000, '\n');
            cin>>choice;
            fail = cin.fail();
            if (!fail)
            {
                if (choice >6 || choice <1)
                {
                    fail = true;
                }
            }
        }

        if(choice==1){
            string src = "Marston";
            string destination ="Potato";
            cout<<src;

            dijkstra(graph, maphash.at(src), maphash.at(destination), hashmap);

            /*cout<<"Enter source name";
            cin>>src;
            cout<<"Enter destination";
            cin>>destination;*/


        }

        else if(choice==2){
            cout << fixed << showpoint;
            cout << setprecision(5);
            double lon;
            double lat;
            string newName;
            int neighborChoice;
            string neighName;
            int neighIndex;
            cout << "\nEnter Longtitude.\n";
            cin >> lon;
            while(cin.fail()){
                cout << "Enter a longtitude number.\n";
                cin >>lon;
            }
            cout << "Enter Latitude.\n";
            cin >> lat;
            while(cin.fail()){
                cout << "Enter a latitude number.\n";
                cin >>lat;
            }
            cout << "Enter location name.\n";
            cin >> newName;
            while(locationNameExists(newName, locations)){
                cout << "Location name already exists. Please enter a different name.\n";
                cin >> newName;
            }
            cout << "Enter the number of neighbors.\n";
            cin >> neighborChoice;
            bool fail = cin.fail();
            if(neighborChoice > locations.size())
                fail = true;
            while(fail){
                cout << "Enter a whole number that is less than the number of locations.\n";
                cin >>neighborChoice;
                if(!cin.fail() && neighborChoice < locations.size())
                    fail = false;
            }
            while(neighborChoice>0){
                cout << "Enter neighbor name.\n";
                cin >> neighName;
                while(!locationNameExists(neighName, locations, neighIndex)){
                    cout << "Neighbor name does not exist, please check spelling and list of neighbors and try again or enter Quit.\n";
                    cin >> neighName;
                    if(neighName=="Quit")
                        break;

                }
                //doesnt word really
                neighborsPassToNode.push_back(locations[neighIndex]);
                distancesPassToNode.push_back(getDist(lon, lat, locations[neighIndex]));
                neighborChoice--;
            }

            locations.push_back(new Node(lon, lat, newName, neighborsPassToNode, distancesPassToNode));
            bool added = false;
            //vector for the adjacency of the new node, which will be pushed to adjMatrix
            vector <int> newLocationAdjacency;


            for (int i = 0; i < locations.size(); ++i){
                for(int j = 0; j< neighborsPassToNode.size(); j++){
                    if (neighborsPassToNode[j]==locations[i])
                    {
                        //pushes 1 to both the vector of the neighbor in adjMatrix and to newLocationAdjaceny if the Nodes match
                        adjMatrix[i].push_back(1);
                        newLocationAdjacency.push_back(1);
                        added = true;
                    }
                }
                //pushes -1 if its not a neighbor
                if(!added && i!=locations.size()-1){
                    adjMatrix[i].push_back(-1);
                    newLocationAdjacency.push_back(-1);
                }
                added = false;
            }

            //pushes 0 for self to the vector for the new node
            newLocationAdjacency.push_back(0);
            //pushes the int vector for the new node onto adjMatrix
            adjMatrix.push_back(newLocationAdjacency);


            neighborsPassToNode.clear();
            distancesPassToNode.clear();
            cout<<"\nLocation Added\n\n";



        }
        else if(choice==3){
            for (int i = 0; i < locations.size(); ++i)
            {
                cout<< "\nName: "<< locations[i]->getName() << " Longtitude: " <<locations[i]->getLon() << " Latitude: "<< locations[i]->getLat();
            }
        }
        else if(choice == 5){
            cout << "\nLocations:\n";
            for (int i = 0; i < locations.size(); ++i)
            {
                cout << " -> " << locations[i]->getName();
                for(int j = 0; j< locations[i]->getNeighbors().size(); j++){
                    cout << "\n    ->  Neighbor: " << locations[i]->getNeighbors()[j]->getName();
                    cout << fixed << showpoint;
                    cout << setprecision(5);
                    cout << "\n    ->  Distance: " << locations[i]->getDistances()[j];
                }
                cout << "\n\n";
            }
        }
        else if(choice == 6){
            return 0;
        }
    }

    return 0;
}
