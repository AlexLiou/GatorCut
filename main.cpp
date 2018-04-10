#include<iostream>
#include<cmath>
#include <string>
#include <fstream>
#include<vector>
#include <sstream>
#include "termproj.h"
using namespace std;



vector<string> createTokens(string s){
//    std::cout<<s<<endl;
    vector <string> lineVector;
    string temp="";
    for (int i=0; i<s.length();i++){
        char c = s[i];
//        std::cout<<c<<" "<<endl;
//        if(c==' '){
//            continue;
//        }
        if(c != ','){
            temp+=c;
//            std::cout<<temp<<endl;
        }
        else{
//            if(temp != " " && !temp.empty()){
//                std::cout<<"THISIS"<<temp<<endl;
                lineVector.push_back(temp);
//            }
            //overload operator
            string jk;
            stringstream ss;
            ss << c;
            ss >> jk;
            lineVector.push_back(jk);
            temp = "";
        }
    }
    return lineVector;
}



int main(){
    vector <Node> Locations;//vector of node locations
    vector <string> temp;
    vector <string> token;
    string filename;
    std::cout << "Please enter the name of the input file: " << endl;
    std::cin >> filename;

    ifstream infile(filename);
//    infile.open(filename);
    string str;
    while(getline(infile,str)){
        temp = createTokens(str);

        for(int i=0; i<temp.size();i++){
            if (temp[i] != ","){
                token.push_back(temp[i]);
            }
    }
    }

    for(auto p: token){
        std::cout<<p<<endl;
    }

    Node* a = new Node("Marston",0.0,0.0);
    Node* b = new Node("Turlington",2.0,3.0);
    double nodeDist= getDist(*a,*b);
    std::cout<<nodeDist<<endl;
    return 0;
}

