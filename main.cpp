#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <cmath>
#include <stdlib.h> //needed for atof -> convert string to double
#include "main.h"
using namespace std;

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

LinkedNode::LinkedNode(double distance, Node* neighbor){
this->distance = distance;
this->neighbor = neighbor;
this->next=NULL;
}
double getDist(double lat1d, double lon1d, double lat2d, double lon2d){

    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r)/2);
    v = sin((lon2r - lon1r)/2);
    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
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
    for (int i=0; i<s.length();i++)
    {
        char c = s[i];
//      std::cout<<c<<" "<<endl;
//        if(c==' '){
//            continue;
//        }
        if(c != ',')
        {
            temp+=c;
           //std::cout<<temp<<endl;
        }
        else
        {
//          if(temp != " " && !temp.empty())
//          {
//                std::cout<<"THISIS"<<temp<<endl;
                lineVector.push_back(temp);
//          }
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
bool locationNameExists(string checkCase, vector <Node*> locations, int & index)
{
    for(int i=0; i<locations.size(); i++)
    {
        if(checkCase==locations[i]->getName())
        {
            index = i;
            return true;
        }
    }
    return false;
}
//checks to see if a string matches the name of a node in the locations vector
bool locationNameExists(string checkCase, vector <Node*> locations)
{
    for(int i=0; i<locations.size(); i++)
    {
        if(checkCase==locations[i]->getName())
            return true;
    }
    return false;
}

int main()
{
    vector <Node*> locations;//vector of node locations
    vector <string> temp;
    vector <string> token;
    vector <vector <int> > adjMatrix;
    vector <LinkedNode*> adjList;
    vector <Node*> neighborsPassToNode;
    vector <double> distancesPassToNode;
    int lines = 0;
    string filename = "map.txt";

    ifstream infile(filename);
	//    infile.open(filename);
    if (!infile)
    {
        printf("Failed to open map.txt");
        return 0;
    }
    string str;
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
    string name;
    double lon;
    double lat;
    int iter = 0;	
    //each line correlates to a new node so each loop creates a node
    for(int i = 1; i <= lines; i++)
    {
    	name = token[iter];
    	iter++;
    	lon = atof(token[iter].c_str()); //converts string to double
    	iter++;
    	lat= atof(token[iter].c_str());
    	if(i!=lines)
        {
    		iter++;
        }
    	locations.push_back(new Node(lat, lon, name, neighborsPassToNode, distancesPassToNode));
    }
    //prints out to test the Node locations vector
    /*for(int i= 0; i<locations.size(); i++){
    	cout<< locations[i]->getlon();
    	cout<< locations[i]->getY();
    	cout<< locations[i]->getName();
    }*/

    string adjacencyFileName = "adj.txt";
    //clears token vector for reuse
    token.clear();
    ifstream afile(adjacencyFileName);
    if (!afile)
    {
        printf("Failed to open adj.txt");
        return 0;
    }
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
for(int i = 0; i < locations.size(); i++){
        for (int j = 0; j < locations[i]->getNeighbors().size(); j++)
        {
            //if it is the first neighbor, create the head of the linked list in the adjList
           if(j==0){
                adjList.push_back(new LinkedNode(locations[i]->getDistances()[j], locations[i]->getNeighbors()[j]));
                tempLink = adjList[i];
                //cout << "\n Pushed " << locations[i]->getNeighbors()[j]->getName() << " to adjList at the index for " << locations[i]->getName();
            }
            //if adjList has already been pushed the head at this index, set the next node to as the j neighbor
            else{
                tempLink->setNext(new LinkedNode(locations[i]->getDistances()[j], locations[i]->getNeighbors()[j]));
                tempLink=tempLink->getNext();
                //cout << "\n Added " << locations[i]->getNeighbors()[j]->getName() << " to adjList at the index for " << locations[i]->getName();
            }
        }
}

//tests adjList
/*
for(int i = 0; i < adjList.size(); i++){
    tempLink = adjList[i];
        for (int j = 0; j < locations[i]->getNeighbors().size(); ++j)
        {
            cout << locations[i]->getName() << endl;
            cout << "   Neighbor: " << tempLink->getNeighbor()->getName() << endl; 
            cout << "       Distance: " << tempLink->getDistance() << endl;
            tempLink = tempLink->getNext();
            
        }
}
*/
    int V = lines;
    struct Graph* graph = createGraph(V);
    for(int i = 0; i < adjList.size(); i++){
    tempLink = adjList[i];
        for (int j = 0; j < locations[i]->getNeighbors().size(); ++j)
        {
            cout << locations[i]->getName() << endl;
            cout << "   Neighbor: " << tempLink->getNeighbor()->getName() << endl; 
            cout << "       Distance: " << tempLink->getDistance() << endl;
            addEdge(graph, )
            tempLink = tempLink->getNext();

            
        }
    }

//prints neighbors
/*for(int i=0;i<lines;i++)  
    {
        for(int j=0;j<locations[i]->getNeighbors().size();j++)  
        {
            cout<<locations[i]->getNeighbors()[j]->getName(); 
        } 
    }
*/



    int choice;
    while(choice!=6)
    {
        cout<<"\n--------------------Choices--------------------\n";
        cout<<"                 1. Directions\n";
        cout<<"                2. Add Location\n";
        cout<<"               3. Print Locations\n";
        cout<<"                5. Database test\n";
        cout<<"                    6. Exit\n\n";
        cin >> choice;
        bool fail = cin.fail();
        if(choice > 6 || choice < 1)
        {
            fail = true;
            //validates the user input for being both an integer and in the range of valid input
        }
        while(fail)
        {
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

        if(choice==1)
        {
            cout<<"Not yet implemented.\n";
        } else if(choice==2){
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
        locations.push_back(new Node(lon, lat, newName, neighborsPassToNode, distancesPassToNode));
        cout << "Enter the number of neighbors.\n";
        cin >> neighborChoice;
        bool fail = cin.fail();
        if(neighborChoice > locations.size())
            fail = true;
        while(fail){
            cout << "Enter a whole number that is less than half the number of locations.\n";
            cin >>neighborChoice;
            if(!cin.fail() && neighborChoice < locations.size()/2)
                fail = false;
        }
        //takes in a neighbor name and adds it to adjList and locations
        while(neighborChoice>0){
            cout << "Enter neighbor name.\n";
            cin >> neighName;
            while(!locationNameExists(neighName, locations, neighIndex)){
                cout << "Neighbor name does not exist, please check spelling and list of neighbors and try again or enter Quit.\n";
                cin >> neighName;    
                if(neighName=="Quit")
                    break;

            }
           double dist;
        for(int i=0; i<locations.size()-1; i++){
     
            if(neighName == locations[i]->getName()){
                //locations update
                dist = getDist(locations[i]->getLat(), locations[i]->getLon(), locations[locations.size()-1]->getLat(), locations[locations.size()-1]->getLon());
                locations[i]->pushNeighbors(locations[locations.size()-1]);
                locations[i]->pushDistances(dist);
                locations[locations.size()-1]->pushNeighbors(locations[i]);
                locations[locations.size()-1]->pushDistances(dist);
                
                //adjList update
                
                //updates the neighbor -> iterates through the list and links the end to the new neighbor
                tempLink = adjList[i];
                for(int j = 0; j < locations[i]->getNeighbors().size()-2; j++){
                    tempLink=tempLink->getNext();
                }
                tempLink->setNext(new LinkedNode(dist, locations[locations.size()-1]));
                
                //updates the new Node -> pushes to adjList if an element of a adjList hasnt been made for this index
                if(adjList.size()<locations.size()){
                    adjList.push_back(new LinkedNode(dist, locations[i]));
                }
                //otherwise iterates through and adds a linked node for the new neighbor
                else{
                    tempLink=adjList[adjList.size()-1];
                    while(tempLink->getNext()!=NULL){
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
        
        
        cout<<"\nLocation Added\n\n";



    }else if(choice==3)
        {
            for (int i = 0; i < locations.size(); ++i)
            {
                cout<< "\nName: "<< locations[i]->getName() << " Longtitude: " <<locations[i]->getLon() << " Latitude: "<< locations[i]->getLat();
            }
        } else if(choice == 5)
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
                    cout << "\n    ->  Distance: " << locations[i]->getDistances()[j];
                }
                cout << "\n\n";
            }
        } else if(choice == 6)
        {
        return 0;
        }
    }
    return 0;
}
