#include<iostream>
#include<cmath>
#include <string>
#include <fstream>
#include<vector>
#include <sstream>
#include <stdlib.h> //needed for atof -> convert string to double
#include "main.h"
using namespace std;

Node::Node(int x, int y, string name){
	this->x=x;
	this->y=y;
	this->name = name;
}

string Node::getName(){
	return this->name;
}

double Node::getX(){
	return this->x;
}

double Node::getY(){
	return this->y;
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



int main(){
    vector <Node*> locations;//vector of node locations
    vector <string> temp;
    vector <string> token;
    vector <vector <int> > adjMatrix;
    int lines = 0;
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
    }*/
    //iterates through the token list and creates a Node for each one specified
    string name;
    double x;
    double y;
    int iter = 0;	
    //each line correlates to a new node so each loop creates a node
    for(int i = 1; i <= lines; i++){
    	name = token[iter];
    	iter++;
    	x = atof(token[iter].c_str()); //converts string to double
    	iter++;
    	y= atof(token[iter].c_str());
    	if(i!=lines)
    		iter++;
    	locations.push_back(new Node(x, y, name));
    }
    //prints out to test the Node locations vector
    /*for(int i= 0; i<locations.size(); i++){
    	cout<< locations[i]->getX();
    	cout<< locations[i]->getY();
    	cout<< locations[i]->getName();
    }*/

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
    }*/


    /*Node* a = new Node("Marston",0.0,0.0);
    Node* b = new Node("Turlington",2.0,3.0);
    double nodeDist= getDist(*a,*b);
    std::cout<<nodeDist<<endl;
    */
    return 0;
}









