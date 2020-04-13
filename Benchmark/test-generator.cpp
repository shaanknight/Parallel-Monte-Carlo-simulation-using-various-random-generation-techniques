#include<iostream>
#include <random>
#include <vector>
#include <stack>
#include "testlib.h"
using namespace std;

constexpr double PI = 3.14159265358979323846264338;
typedef unsigned long long ll;

// Struct to hold positions
struct Point {
    double x,y;
    Point(){}
    Point(double x, double y){
        this->x = x;
        this->y = y;
    }
};

// A global point needed for sorting points with reference to the first point
Point base;

// A utility function to find next to top in a stack
Point nextToTop(stack<Point> &S){
	Point p = S.top();
	S.pop();
	Point res = S.top();
	S.push(p);
	return res;
}

double distSq(Point p1, Point p2){
	return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}

int orientation(Point p, Point q, Point r){
	double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if(abs(val) < 0.0000000000005) return 0;
	return (val > 0)? 1: 2; // clock or counterclock wise
}

// A function used by library function qsort() to sort an array of
// points with respect to the first point
int compare(const void *vp1, const void *vp2){
    Point *p1 = (Point *)vp1;
    Point *p2 = (Point *)vp2;

    // Find orientation
    int o = orientation(base, *p1, *p2);
    if (o == 0)
    	return (distSq(base, *p2) >= distSq(base, *p1))? -1 : 1;

    return (o == 2)? -1: 1;
}

// Prints convex hull of a set of n points.
std::vector<Point> convexHull(Point points[], int n){
    // Find the bottommost point
    double ymin = points[0].y;
    int min = 0;
    for (int i = 1; i < n; i++)
    {
    	double y = points[i].y;

    	// Pick the bottom-most or chose the left
    	// most point in case of tie
    	if ((y < ymin) || (ymin == y &&
    		points[i].x < points[min].x))
    		ymin = points[i].y, min = i;
    }

    // Place the bottom-most point at first position
    base = points[min];
	points[min] = points[0];
	points[0] = base;
    qsort(&points[1], n-1, sizeof(Point), compare);

    int m = 1; // Initialize size of modified array
    for (int i=1; i<n; i++)
    {
    	while (i < n-1 && orientation(base, points[i],points[i+1]) == 0)
    		i++;
    	points[m] = points[i];
    	m++; // Update size of modified array
    }

    std::vector<Point> vertices;
    if (m < 3) return vertices;

    stack<Point> S;
    S.push(points[0]);
    S.push(points[1]);
    S.push(points[2]);
    for (int i = 3; i < m; i++)
    {
    	while (orientation(nextToTop(S), S.top(), points[i]) != 2)
            S.pop();
    	S.push(points[i]);
    }
    while (!S.empty()){
    	Point p = S.top();
        vertices.push_back(p);
    	S.pop();
    }
    return vertices;
}

int main(int argc, char* argv[])
{
    registerGen(argc, argv, 1);

    const int mod = 100000.0;
    Point points[1024];
    int n = 40;
    for(int i=0;i<n;i++){
        points[i] = Point((double)(rnd.next(0,mod-1))/(mod-1.0),(double)(rnd.next(0,mod-1))/(mod-1.0));
    }
    std::vector<Point> vertices = convexHull(points,n);
    for(int i=0;i<(int)vertices.size();i++)
        cout << vertices[i].x << " " << vertices[i].y << "\n";
    return 0;
}