#include <string>
using namespace std;
class Node{

private:
	double x;
	double y;
	string name;

public:
	string getName();
	double getX();
	double getY();
	Node(int x, int y, string name);
};


/*void PrintMenu(*Node [] locations);
void addNodeLocation(string name, int x, int y);
string [] Dijkstra(*Node [] locations);
*/