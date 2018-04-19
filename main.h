#include <string>
#include <math.h>
#include <cmath>
#define earthRadiusKm 6371.0
using namespace std;



class Node{

private:
	double lon;
	double lat;
	string name;
	vector <Node*> neighbors;
	vector <double> distances;

public:
	string getName();
	double getLon();
	double getLat();
	vector <Node*> getNeighbors();
	vector <double> getDistances();
	void pushNeighbors(Node* neighbor);
	void pushDistances(double distance);
	Node(double lon, double lat, string name, vector <Node*> neighbors, vector <double> distances);
};

double deg2rad(double deg) {
  return (deg * M_PI / 180);
}

class LinkedNode{

		private:
	double distance;
	Node* neighbor;
	LinkedNode* next;
		public:
	double getDistance();
	double getDist(double lat1d, double lon1d, double lat2d, double lon2d);
	Node* getNeighbor();
	LinkedNode* getNext();
	void setNext(LinkedNode* next);
	LinkedNode(double distance, Node* neighbor);
};


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

double getDist(double lonOne, double latOne, Node* secondNode){
	double lat1r, lon1r, lat2r, lon2r, u, v;
	lat1r = latOne * M_PI/180;
	lon1r = lonOne * M_PI/180;
	lat2r = secondNode->getLat() * M_PI/180;
	lon2r = secondNode->getLon() * M_PI/180;
  	u = sin((lat2r - lat1r)/2);
  	v = sin((lon2r - lon1r)/2);
  return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));

}
