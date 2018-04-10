//
// Created by Asia Clifford on 4/10/18.
//
#include <iostream>
#include <cmath>
using namespace std;
#ifndef PROG2TERMPROJ_TERMPROJ_H
#define PROG2TERMPROJ_TERMPROJ_H

#endif //PROG2TERMPROJ_TERMPROJ_H

class Node{
private:
    string name;
    double xval,yval;
public:
    //constructor
    Node(string n,double x=0.0,double y=0.0){
        name =n;
        xval=x;
        yval=y;
    }
    //get methods
    string getName(){return name;}
    double getX(){return xval;}
    double getY(){return yval;}
};
double getDist(Node c,Node d){
    double x1= c.getX();
    double x2=d.getX();
    double y1=c.getY();
    double y2=d.getY();

    double xDiff = x1-x2;
    double yDiff = y1-y2;
    double dist= sqrt(pow(xDiff,2)+pow(yDiff,2));
    return dist;
}