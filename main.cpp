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

// Constructor for LinkedNode Class:
// LinkedNode Class, used for data representation in LinkedList.
// Sets Distance and Neighbor equal to parameter.
// Sets next equal to NULL.
LinkedNode::LinkedNode(double distance, Node* neighbor)
{
    this->distance = distance;
    this->neighbor = neighbor;
    this->next = NULL;
}

// Utility Function for LinkedNode Class:
// Returns doubel distance variable in LinkedNode Class.
double LinkedNode::getDistance()
{
    return this->distance;
}

// Utility Function for LinkedNode Class:
// Returns Node Pointer Neighbor in LinkedNode Class
Node* LinkedNode::getNeighbor()
{
    return this->neighbor;
}

// Utility Function for LinkedNode Class:
// Returns Next Node Pointer in LinkedNode Class
LinkedNode* LinkedNode::getNext()
{
    return this->next;
}

// Utility Function for LinkedNode Class:
// Sets Next Node Pointer in LinkedNode Class
void LinkedNode::setNext(LinkedNode* next)
{
    this->next = next;
}

// Constructor for Node Class:
// Used for data representation in Vector of Nodes for File Input. 
// Sets Longitude, Latitude, name, neighbors, and distances equal to the parameters. 
// 
Node::Node(double lon, double lat, string name, vector <Node*> neighbors, vector <double> distances)
{
    this->lon=lon;
    this->lat=lat;
    this->name = name;
    this->neighbors=neighbors;
    this->distances=distances;
}

// Utility Function for Node Class:
// Returns string name variable in Node Class
string Node::getName(){
    return this->name;
}

// Utility Function for Node Class:
// Returns double longitude variable in Node Class
double Node::getLon()
{
    return this->lon;
}

// Utility Function for Node Class:
// Returns double latitigude variable in Node Class
double Node::getLat()
{
    return this->lat;
}

// Utility Function for Node Class:
// Returns Vector of Node Pointers in Node Class.
vector <Node*> Node::getNeighbors()
{
    return this->neighbors;
}

// Utility Function for Node Class:
// Pushes Node* to end of Vector of Node Pointers representing Neighbors in Node Class.
void Node::pushNeighbors(Node* neighbor)
{
    this->neighbors.push_back(neighbor);
}

// Utility Function for Node Class:
// Returns Vector of doubles representing distances of Neighbors in Node Class.
vector <double> Node::getDistances()
{
    return this->distances;
}

// Utility Function for Node Class:
// Pushes double distance to end of Vector of doubles representing distances of Neighbors in Node Class.
void Node::pushDistances(double distance)
{
    this->distances.push_back(distance);
}

// File I/O Function for Main:
// Parses Tokens from File Input, separates
// File lines by "," delimiter.
vector<string> createTokens(string s)
{
//    std::cout<<s<<endl;
    vector <string> lineVector;
    string temp="";
    // Loops through string s, checking if each character in the string is the delimiter. 
    // If the character is not the delimiter, adds the character to the buffer string "temp".
    // If the character is the delimiter, pushes the buffer string to the lineVector Vector and erases buffer string.
    for (int i=0; i<s.length();i++)
    {
        // Character to be inspected.
        char c = s[i];
        //       std::cout<<c<<" "<<endl;
//        if(c==' ')
//            continue;
//      
        // Checks if character is delimiter.
        if(c != ',')
        {
            temp+=c;
            //std::cout<<temp<<endl;
        }
        else
        {
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

// Utility Function for Main class:
// Checks to see if a string matches the name of a node in the locations vector.
// Changes the value of index to the index where a match is found
bool locationNameExists(string checkCase, vector <Node*> locations, int & index)
{
    for(int i=0; i<locations.size(); i++)
    {
        if(checkCase==locations[i]->getName())
        {
            index = i;
        } 
        return true;
    }
    return false;
}

// Utility Function for Main class:
// Checks to see if a string matches the name of a node in the locations vector
bool locationNameExists(string checkCase, vector <Node*> locations)
{
    for(int i=0; i<locations.size(); i++)
    {
        if(checkCase==locations[i]->getName())
        {
                return true; 
        }
    }
    return false;
}

// Structure for Adjacency List Node.
// Represents Adjacency List Node with int destination for index, double weight for 
// distance between Node and Destination, and Node Pointer for connected Node.
struct AdjListNode
{
    int dest;
    double weight;
    struct AdjListNode* next;
};

// Structure for Adjacency List
// Represents an Array of LinkedLists, with a head Node Pointer AdjListNode.
struct AdjList
{
    struct AdjListNode *head;  // pointer to head node of list
};

// Structure for Graph
// Represents an Array of Adjacency Lists.
// Size of array will be V (number of vertices in graph)
struct Graph
{
    int V;
    struct AdjList* array;
};

// Utility Function for Adjacency List Node Class:
// Creates a new Adjacency List Node, allocating memory for size of AdjListNode.
// Sets new nodes destination, weight and next pointer equal to parameters.
struct AdjListNode* newAdjListNode(int dest, double weight)
{
    struct AdjListNode* newNode =
            (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}

// Utility Function for Graph Structure:
// Creates Graph, allocating memeory for it, based on number of Vertexes V. 
struct Graph* createGraph(int V)
{
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;

    // Create an array of adjacency lists.  Size of array will be V
    graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));

    // Initialize each adjacency list as empty by making head as NULL
    for (int i = 0; i < V; ++i){
             graph->array[i].head = NULL; 
    }

    return graph;
}

// Utility Function for Graph Structure:
// Addas an Edge, an Adjacency List Node to the Graph's array of Adjacency Lists.
// Adjacency List Node represents the physical distance between two Nodes within the Graph.
// Ex: The distance between two buildings is represented as an Edge.
void addEdge(struct Graph* graph, int src, int dest, double weight)
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

// Constructur for MinHeapNode Structure.
// Declares parameters v and dist, used to represent vector and distance.
struct MinHeapNode
{
    int  v;
    double dist;
};

// Constructor for MinHeap Structure
// Declares size, capacity, and position pointer.
struct MinHeap
{
    int size;      // Number of heap nodes present currently
    int capacity;  // Capacity of min heap
    int *pos;     // This is needed for decreaseKey()
    struct MinHeapNode **array;
};

// Utility Function for MinHeapNode Structure
// Creates a new Min Heap Node, allocating memeory for it and initializing the v and dist.
struct MinHeapNode* newMinHeapNode(int v, double dist)
{
    struct MinHeapNode* minHeapNode =
            (struct MinHeapNode*) malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}

// Utility Function for MinHeap Structure
// Creates the MinHeap, allocating memory for it and initializing the position, size,
// capacity, and array for MinHeapNodes.
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

// Utility Function for MinHeapNodes
// Swaps two nodes with each other in each position. Used in the Heapify Function.
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b)
{
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// Utility Function for MinHeap:
// A standard function to heapify at given index/idx. Keeps the smallest MinHeapNode on top of the MinHeap.
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

// Utility Function for MinHeap:
// Checks if the given minHeap is empty or not. Returns int value. 
int isEmpty(struct MinHeap* minHeap)
{
    return minHeap->size == 0;
}

// Utility Function for MinHeap:
// Standard function to extract minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap* minHeap)
{
    // Returns NULL if nothing is in MinHeap.
    if (isEmpty(minHeap))
    {
        return NULL;
    }

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

// Utility Function for MinHeap:
// Function decreases dist value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap* minHeap, int v, double dist)
{
    // Gets index of v in MinHeap Array.
    int i = minHeap->pos[v];

    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;

    // Travel up while the complete tree is not heapified.
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

// Utility Function for MinHeap:
// Checks if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct MinHeap *minHeap, int v)
{
    if (minHeap->pos[v] < minHeap->size)
    {
        return true;
    }
    return false;
}

// Utility Function for dijstra's Algorithm:
// Recursive Function that prints the nodes in order from source to destination.
void printPath(int parent[], int src, int destination, map<int, string> hashmap)
{
    // Base Case
    if(destination==src) return;
  
    // Recursive Check
    printPath(parent, src, parent[destination], hashmap);
    cout<<hashmap.at(destination)<<"\t\t";
    destination=parent[destination];
}

// Utility Function for Dijstra's Algorithm:
// Used to print the array of answers for the Dijstra's Algorithm. Gives source, destination, 
// distance in kilometers and steps. Calls Recursive Function printPath to list order of nodes
// from source to destination. Uses hashmap to convert indexes obtained from Dijkstra's Algorithm
// and converts them into building names.
void printArr(double dist[], int src, int parent[], int destination, map<int, string> hashmap, string transp)
{
    
    printf("Source:\t\tDestination:\t\tDistance from Source:\n");

    cout<<hashmap.at(src)<<"\t\t"<<hashmap.at(destination)<<"\t\t"<<dist[destination]<< " km";

        cout << fixed << showpoint;
        cout << setprecision(3);

    if(transp == "Walk"){
        double steps = dist[destination]*1320;
        cout<<"\t=\t"<<steps<<" steps";
        cout<<endl;
        double time = dist[destination]/0.1166667;
        cout<<endl;
        cout<<"Time to destination (walking): "<<time<< " minute(s)";
    }

    if(transp =="Bike"){
        cout<<endl;
        double time = dist[destination]/0.3333333;
        cout<<endl;
        cout<<"Time to destination (biking): "<<time<<" minute(s)";
    }

    if(transp =="Skateboard"){
        cout<<endl;
        double time = dist[destination]/0.2;
        cout<<endl;
        cout<<"Time to destination (skateboarding): "<<time<<" minute(s)";
    }

    if(transp=="Long-board"){
        cout<<endl;
        double time = dist[destination]/0.25;
        cout<<endl;
        cout<<"Time to destination (longboarding): "<<time<<" minute(s)";
    }


    cout << endl;
    cout << "\nPath from " << hashmap.at(src) << endl;
    printPath(parent, src, destination, hashmap);

}

// Dijstra's Algorithm: 
// The main function that calulates distances of shortest paths from source node to destination node
// vertices. It is a O(ELogV) function. Calls the printArr function to print list of Nodes. 
double dijkstra(struct Graph* graph, int src,int dest, map<int, string> hashmap, string transp, bool food)
{   
  
  //starting
    if(food == true){ 
        int V = graph->V;// Get the number of vertices in graph
    double dist[V];      // dist values used to pick minimum weight edge in cut
    int destinIndex= dest;
    int parent [V];
      
    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = std::numeric_limits<double>::max();
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
            if (isInMinHeap(minHeap, v) && dist[u] != std::numeric_limits<double>::max() &&
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
         double x = dist[dest];
        return x;
    }else{
    int V = graph->V;// Get the number of vertices in graph
    double dist[V];      // dist values used to pick minimum weight edge in cut
    int destinIndex= dest;
    int parent [V];
    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = std::numeric_limits<double>::max();
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
            if (isInMinHeap(minHeap, v) && dist[u] != std::numeric_limits<double>::max() &&
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
    printArr(dist, src, parent, dest, hashmap, transp);
    return 0;
    }
}


int main()
{
    //Vector of node locations: Holds data to be printed in methods to check database.
    vector <Node*> locations;
  
    //Temporary Vector.
    vector <string> temp;
  
    //Vector of Tokens, taken from file string.
    vector <string> token;
  
    //used to figure out which nodes are neighbors
    vector <vector <int> > adjMatrix; 
    vector <LinkedNode*> adjList;
    vector <Node*> neighborsPassToNode;
    vector <double> distancesPassToNode;
    int lines = 0;
  
    // Calls name of File.
    string filename = "map.txt";
        // Input File of Building names and longitude and latitude coordinates.
    // Data is separated with commas.
    ifstream infile(filename);
    //    infile.open(filename);
    string str;
  
    // Checks if file is inputted successfully.
    if (!infile)
    {
      cout << "Error in accessing file.";
      return 0;
    }
  
    // Loops through map.txt, creating vector of tokens.
    while(getline(infile,str))
    {
        temp = createTokens(str);
        lines++;
        for(int i=0; i<temp.size();i++)
        {
            //if (temp[i] != ","){
                token.push_back(temp[i]);
           // }
            }
    }
  
    //tests token vector
    /*for(auto p: token){
        std::cout<<p<<endl;
    }*/
    //iterates through the token list and creates a Node for each one specified
  
    // Variables to represent data in tokens.
    string name;
    double lon;
    double lat;
    int iter = 0;
  
    //each line correlates to a new node so each loop creates a node
    for(int i = 1; i <= lines; i++)
    {
        name = token[iter];
        iter++;
      
        //converts string to double
        lon = atof(token[iter].c_str()); 
        iter++;
      
        //converts string to double
        lat= atof(token[iter].c_str());
      
        if(i!=lines)
        {
            iter++;
        }
        // Adds new Node to the Vector of locations, represeting the building.
        locations.push_back(new Node(lat, lon, name, neighborsPassToNode, distancesPassToNode));
    }
  
    //prints out to test the Node locations vector
    /*for(int i= 0; i<locations.size(); i++){
        cout<< locations[i]->getlon();
        cout<< locations[i]->getY();
        cout<< locations[i]->getName();
    }*/

    // File Name for adjacency matrix for map.txt.
    // 0 represents building is itself.
    // 1 represents building is a neighbor.
    // -1 represents building is not a neighbor.
  
    string adjacencyFileName = "adj.txt";
  
    //clears token vector for reuse
    token.clear();
    ifstream afile(adjacencyFileName);
  
    // Loops through file, getting tokens based on the comma delimiter.
    while(getline(afile,str))
    {
        temp = createTokens(str);
        for(int i=0; i<temp.size();i++)
        {
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
    for(int i = 1; i <= lines; i++)
    {
        for(int j = 0; j<lines; j++)
        {
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
    }*/

    // Creates the Vector of Neighbors and corresponding Distances for each location Node. 
    for(int i = 0; i < lines; i++)
    {
        for(int j = 0; j<lines; j++)
        {
            if(adjMatrix[i][j]==1)
            {
                locations[i]->pushNeighbors(locations[j]);
                locations[i]->pushDistances(getDist(locations[i]->getLat(), locations[i]->getLon(), locations[j]->getLat(), locations[j]->getLon()));
            }

        }
    }

    //creates the adjacency list by making a node for each neighbor at a locations index
    LinkedNode* emptyLink;
    LinkedNode* tempLink;
  
    for(int i = 0; i < locations.size(); i++)
    {
        for (int j = 0; j < locations[i]->getNeighbors().size(); j++)
        {
        //if it is the first neighbor, create the head of the linked list in the adjList
                if(j==0)
            {
                adjList.push_back(new LinkedNode(locations[i]->getDistances()[j], locations[i]->getNeighbors()[j]));
                tempLink = adjList[i];
                //cout << "\n Pushed " << locations[i]->getNeighbors()[j]->getName() << " to adjList at the index for " << locations[i]->getName();
            }
            //if adjList has already been pushed the head at this index, set the next node to as the j neighbor
            else
            {
                tempLink->setNext(new LinkedNode(locations[i]->getDistances()[j], locations[i]->getNeighbors()[j]));
                tempLink=tempLink->getNext();
                //cout << "\n Added " << locations[i]->getNeighbors()[j]->getName() << " to adjList at the index for " << locations[i]->getName();
            }
        }
    }

    //tests adjList
    /*for(int i = 0; i < adjList.size(); i++){
        tempLink = adjList[i];
        cout << locations[i]->getName() << endl;
        cout << "   Neighbor: " << tempLink->getNeighbor()->getName() << endl;
        cout << "       Distance: " << tempLink->getDistance() << endl;
            while(tempLink->getNext()!=NULL)
            {
                tempLink = tempLink->getNext();
                cout << locations[i]->getName() << endl;
                cout << "   Neighbor: " << tempLink->getNeighbor()->getName() << endl;
                cout << "       Distance: " << tempLink->getDistance() << endl;
            }
    }*/

    //prints neighbors
    /*for(int i=0;i<lines;i++)
        {
            for(int j=0;j<locations[i]->getNeighbors().size();j++)
            {
                cout<<locations[i]->getNeighbors()[j]->getName();
            }
        }
    */
  
    // Variable V represents the number of Vectors in the locations Vector.
    int V = locations.size();
  
    // Creates the Hashmap for int keys and string values, used for Dijstra's Algorithm
    // to return the related strings for the indexes it finds for the shortest path.
    map<int,string> hashmap;
  
    for (int i=0; i < V; i++)
    {
        hashmap[i] = locations[i]->getName();
    }

    // Creates the Hashmap for string keys and int values, used for user input to 
    // find the corresponding index for string to be used in Dijstra's Algorithm.
    map<string, int> maphash;
  
    for (int i=0; i < V; i++)
    {
        maphash[locations[i]->getName()] = i;
    }
  
    // Creates Graph Structure. Uses addEdge Function to creates Edges representing
    // the relationships between locations and their neighbors.
    struct Graph* graph = createGraph(V);
  
    for (int i=0; i < V; i++)
    {

        for(int j = 0; j< locations[i]->getNeighbors().size(); j++)
        {

            int num = maphash.at(locations[i]->getNeighbors()[j]->getName());
            addEdge(graph, i, num, locations[i]->getDistances()[j]);
            addEdge(graph, num, i, locations[i]->getDistances()[j]);
        }
    }


    // While Loop to continue outputting menu until 6 is entered and loop is ended.
    int choice;
  
    while(choice!=5)
    {
        // Prints menu.
        cout<<endl;
            cout<<"\n--------------------Choices--------------------\n";
            cout<<"                 1. Directions\n";
            cout<<"                2. Add Location\n";
            cout<<"               3. Print Locations\n";
            cout<<"                4. View Database \n";
            cout<<"                    5. Exit\n\n";
            cin >> choice;
      
        // Checks if user input is incorrect for its purpose.
            bool fail = cin.fail();
      
            if(choice > 5 || choice < 1){
                fail = true;
            //validates the user input for being both an integer and in the range of valid input
        }
      
        while(fail)
        {
                cout<<"Invalid input, please input a whole number between 1 and 5.\n";
          
                // Flushes input stream to allow for user input.
              cin.clear();
                cin.ignore(10000, '\n');
                cin>>choice;
                fail = cin.fail();
          
              if (!fail)
              {
                  if (choice >5 || choice <1)
                  {
                      fail = true;
                  }
              }
        }
      
      
        // Option 1: Getting Shortest Path between two buildings
            if(choice==1)
        {
            cout << "\nEnter starting location: ";
            string src;
            cin >> src;

            cout << "\nEnter end location: ";
            string destination;
            cin >> destination;
                        
            // Checks if User Input is invalid, not in list of destinations available.
            if (maphash.find(src) == maphash.end())
            {
                cout << "\nLocation(s) not in database.";
                continue;
            } 
                else if (maphash.find(destination) == maphash.end())
            {
                cout << "\nLocation(s) not in database.";
                continue;
            }
            if (maphash.find(src) == maphash.end())
            {
                cout << "\nLocation(s) not in database.";
                continue;
            } 
            else if (maphash.find(destination) == maphash.end())
            {
                cout << "\nLocation(s) not in database.";
                continue;
            }
          
            string transport[4] = {"Walk", "Bike", "Skateboard", "Long-board"};
            string transp;

            bool cont = true;
          
            while(cont){
                    //users choosed their transportation method which will tell them how long
                    // it will take to reach their destination
                cout <<"\nChoose transportation method:\n";
                cout <<"------------Walk------------\n";
                cout <<"------------Bike------------\n";
                cout <<"---------Skateboard---------\n";
                cout <<"---------Long-board---------\n";
                cout<<endl;
                cin >> transp;

                if(transp!=transport[0] && transp!=transport[1] && transp!=transport[2] && transp!=transport[3]){
                    cout<<endl;
                    cout<<"This transportation method is not available.";
                    cout<<endl;
                    cont=true;
                }
                else{
                    cont=false;
                }
            }
          
            //asks the user if they want to stop for food
            //if no, the program will simply call the algorithm for calculating distance and path between the source and destination
            //if yes, program will discover the fastest way to get to your destination, while making a pit stop for food
            //the program does this by calculating the path through every single food place to your destination
                //takes in user input and checks for incorrect input
            string foodYN;
            cout<<"\nDo you want to stop for food(Y/N)";
            cin>> foodYN;
            bool fail = true;
            if(foodYN == "Y" || foodYN == "N")
                fail = false;

            while(fail){
                cin.clear();
                cin.ignore(1000,'\n');
                cout<< "\nEnter Y or N: ";
                cin >> foodYN; 
                if(foodYN != "Y" && foodYN != "N" && foodYN !="y" && foodYN !="n"){
                    fail = false;
                }
            }
            cout << endl;
            if(foodYN == "Y" || foodYN =="y"){
              
                if (locations[maphash.at(src)]->getName().compare("Reitz-Union") == 0 ||
                    locations[maphash.at(src)]->getName().compare("Little Hall") == 0 ||
                    locations[maphash.at(src)]->getName().compare("Student Rec") == 0 ||
                    locations[maphash.at(src)]->getName().compare("Broward Dining") == 0  ||
                    locations[maphash.at(src)]->getName().compare("Hub") == 0         ||
                    locations[maphash.at(src)]->getName().compare("Hub") == 0         ||
                    locations[maphash.at(src)]->getName().compare("Little-Hall") == 0){

                    cout << "\nThere is food at your source location: " << locations[maphash.at(src)]->getName() << ".\n";
                        dijkstra(graph, maphash.at(src), maphash.at(destination), hashmap, transp, false);
                }else{

                //initlizes variables for tracking and recording the fasted distance food place

                double minDist = 0;
                double trackerDist;
                int bestFood;
              

              //iterates through locations to check what locations have food
                for(int i = 0; i<locations.size(); i++){
                    //if statement -> if the name matches the location has food
                        if (locations[i]->getName().compare("Reitz-Union") == 0 ||
                            locations[i]->getName().compare("Little-Hall") == 0 ||
                            locations[i]->getName().compare("Student Rec") == 0 ||
                            locations[i]->getName().compare("Broward Dining") == 0  ||
                            locations[i]->getName().compare("Hub") == 0         ||
                            locations[i]->getName().compare("Hub") == 0         ||
                            locations[i]->getName().compare("Little-Hall") == 0)

                        {
                    //trackerDist calculates the total distance from src to food to destination
                            
                    trackerDist = dijkstra(graph, maphash.at(src), i, hashmap, transp, true) + dijkstra(graph, i, maphash.at(destination), hashmap, transp, true);
                    //if trackerDist is lower than the previous min and minDist hasnt been changed yet
                    if(trackerDist < minDist || minDist == 0){
                        minDist = trackerDist;
                        bestFood = i;


                    }

                }

                }
              
                dijkstra(graph, maphash.at(src), bestFood, hashmap, transp, false);
                    //prints that there is food at the location
                    cout << "\n\nThere is food at " << locations[bestFood]->getName() << "!!\n\n";
                        //if the destination has food, it does not calculate the zero distance between the two
                if(bestFood!=maphash.at(destination))
                dijkstra(graph, bestFood, maphash.at(destination), hashmap, transp, false);

                }
            }
            else{
                 dijkstra(graph, maphash.at(src), maphash.at(destination), hashmap, transp, false);
            }
        }
      
      // Option 2: Add Location based on required information: Longitude, Latitude, Name, Neighbors,
            else if(choice==2)
      {
                // Prepares for precise coordinate input
            cout << fixed << showpoint;
            cout << setprecision(5);
            double lon;
            double lat;
            string newName;
            int neighborChoice;
            string neighName;
            int neighIndex;
        
                // Asks for Longitude, checks if User input is valid for stated purpose.
            cout << "\nEnter Longtitude.\n";
            cin >> lon;
        
            while(cin.fail())
            {
                    cout << "Enter a longtitude number.\n";
                    cin.clear();
                    cin.ignore(1000,'\n');
                    cin >>lon; 
                }
        
                // Asks for Latitude, checks if User input is valid for stated purpose.
            cout << "Enter Latitude.\n";
            cin >> lat;
        
            while(cin.fail())
            {
                cout << "Enter a latitude number.\n";
                cin.clear();
                 cin.ignore(1000,'\n');
                cin >>lat;
            }
        
                // Asks for Location Name, checks if User input is valid for stated purpose.
            cout << "Enter location name.\n";
            cin >> newName;
        
            while(locationNameExists(newName, locations))
            {
                cout << "Location name already exists. Please enter a different name.\n";
                cin.clear();
                cin.ignore(1000,'\n');
                cin >> newName;
            }
        
            locations.push_back(new Node(lon, lat, newName, neighborsPassToNode, distancesPassToNode));
        
                // Asks for Number of Neighbors, checks if User input is possible and valid for stated purpose.
            cout << "Enter the number of neighbors.\n";
            cin >> neighborChoice;
            bool fail = cin.fail();
        
            if(neighborChoice > locations.size())
            {
                    fail = true; 
            }
        
            while(fail)
            {
                cout << "Enter a whole number that is less than half the number of locations.\n";
                cin >> neighborChoice;
                if(!cin.fail() && neighborChoice < locations.size()/2)
                {
                        fail = false; 
                }
                cin.clear();
                cin.ignore(1000,'\n');
            }
        
            //takes in a neighbor name and adds it to adjList and locations
            while(neighborChoice>0)
            {
                cout << "Enter neighbor name.\n";
                cin >> neighName;
              
                while(!locationNameExists(neighName, locations) || neighName == newName)
                {
                     cin.clear();
                    cin.ignore(1000,'\n');
                    cout << "Neighbor name does not exist, please check spelling and list of neighbors and try again or \nenter Quit if you are done adding neighbors.\n";
                    cin >> neighName;
                    if(neighName=="Quit")
                    {

                        neighborChoice--;
                        break;
                    }
                   
                }
              
            double dist;
              
            // Updates Vector of Neighbors and Distances in Locations based on new Node added.
            for(int i=0; i<locations.size()-1; i++)
            {
                if(neighName == locations[i]->getName())
                {
                    //locations update
                    dist = getDist(locations[i]->getLat(), locations[i]->getLon(), locations[locations.size()-1]->getLat(), locations[locations.size()-1]->getLon());
                    locations[i]->pushNeighbors(locations[locations.size()-1]);
                    locations[i]->pushDistances(dist);
                    locations[locations.size()-1]->pushNeighbors(locations[i]);
                    locations[locations.size()-1]->pushDistances(dist);

                    //adjList update

                    //updates the neighbor -> iterates through the list and links the end to the new neighbor
                    tempLink = adjList[i];
                  
                    for(int j = 0; j < locations[i]->getNeighbors().size()-2; j++)
                    {
                        tempLink=tempLink->getNext();
                    }
                  
                    tempLink->setNext(new LinkedNode(dist, locations[locations.size()-1]));

                    //updates the new Node -> pushes to adjList if an element of a adjList hasnt been made for this index
                    if(adjList.size()<locations.size())
                    {
                        adjList.push_back(new LinkedNode(dist, locations[i]));
                    }
                  
                    //otherwise iterates through and adds a linked node for the new neighbor
                    else
                    {
                        tempLink=adjList[adjList.size()-1];
                        while(tempLink->getNext()!=NULL)
                        {
                            tempLink=tempLink->getNext();
                        }
                        tempLink->setNext(new LinkedNode(dist, locations[i]));
                    }
                }
            }
                neighborChoice--;
            }
        
            //tests adjList
            /*cout << "\n\n";
            for(int i = 0; i < adjList.size(); i++){
                tempLink = adjList[i];
            for (int j = 0; j < locations[i]->getNeighbors().size(); ++j)
            {
                cout << locations[i]->getName() << endl;
                cout << "   Neighbor: " << tempLink->getNeighbor()->getName() << endl;
                cout << "       Distance: " << tempLink->getDistance() << endl;
                tempLink = tempLink->getNext();
            }
    }*/
        
                // Updates Graph, Hashmap, and MapHash for new Node added.
                // Adds to Hashmap and Maphash. Creates Graph anew with updated Locations Vector.
            V++;
            int num = locations.size()-1;
            hashmap[num] = locations[num]->getName();
            maphash[locations[num]->getName()] = num;
            graph = createGraph(V);
        
            for (int i=0; i < V; i++)
            {
                for(int j = 0; j< locations[i]->getNeighbors().size(); j++)
                {
                    int num = maphash.at(locations[i]->getNeighbors()[j]->getName());
                    addEdge(graph, i, num, locations[i]->getDistances()[j]);
                }
            }


            cout<<"\nLocation Added\n\n";


        }
      
      
        // Option 3: Lists the Locations and their coordinate.
            else if(choice==3)
            {
                    for (int i = 0; i < locations.size(); ++i)
                    {
                        cout<< "\nName: "<< locations[i]->getName() << " Longtitude: " <<locations[i]->getLon() << " Latitude: "<< locations[i]->getLat();
                    }
            }
      
      
        // Option 4: Lists the Locations and their Neighbors with corresponding distances in kilometers.
            else if(choice == 4)
        {
                    cout << "\nLocations:\n";
          
                    for (int i = 0; i < locations.size(); ++i)
                    {
                            cout << " -> " << locations[i]->getName();
                
                        for(int j = 0; j< locations[i]->getNeighbors().size(); j++)
                    {
                            cout << "\n    ->  Neighbor: " << locations[i]->getNeighbors()[j]->getName();
                            cout << fixed << showpoint;
                            cout << setprecision(5);
                            cout << "\n    ->  Distance: " << locations[i]->getDistances()[j] << " km";
                        }
                        cout << "\n\n";
                     }
            }
      
            // Option 6: Ends Loop and ends program.
            else if(choice == 5)
          {
              return 0;
                }
    }
  
    return 0;
}
