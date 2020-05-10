#ifndef QUADTREE_HH
#define QUADTREE_HH
#include <vector>

// point structure, contains coordinates and phase
struct Point
{
    double x, y, phase;

    Point(double x = 0, double y = 0, double phase = 0): x(x), y(y), phase(phase) {};
};

// bounding box with center and half side, initializes to square of radius 2 centered at origin
struct Box
{
    Point center;
    double radius;

    Box(Point center = Point(), double radius = 2): center(center), radius(radius) {};

    // checks if p is within bounding box
    bool contains(Point p)
    {
        return p.x<center.x+radius && p.x>center.x-radius && p.y<center.y+radius && p.y>center.y-radius;
    }
};

// each QuadTree contains one point, initializes to square of radius 2 centered at origin
class QuadTree
{
    // four children
    QuadTree* nW;
    QuadTree* nE;
    QuadTree* sW;
    QuadTree* sE;

    // bounding box to represent the boundaries of this quad tree
    Box boundary;

    // centroid, contains coordinate and phase
    // if there are no children, then centroid is the same as the point contained
    // initializes to Point(100,0,0)
    Point centroid;

    // "mass", aka how many points
    int mass;

    // create four children
    void subdivide();

    public:
        // constructor
        QuadTree(Box b);

        // destructor
        ~QuadTree();

        //insert point into quadtree
        bool insert(Point p);

        // update centroid by traversing tree
        void update_centroid();

        std::vector<double> get_centroids(double x, double y, double theta);
};

#endif
