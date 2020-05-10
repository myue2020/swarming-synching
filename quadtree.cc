#include <cmath>
#include <vector>

// point structure, contains coordinates and phase
struct Point {
    double x, y, phase;
    Point(double x = 0, double y = 0, double phase = 0):
        x(x), y(y), phase(phase) {}
};

// bounding box with center and half side, initializes to square of radius 2 centered at origin
struct Box {
    Point center;
    double radius;

    Box(Point center = Point(), double radius = 2):
        center(center), radius(radius) {}

    // checks if p is within bounding box
    bool contains(Point p) {
        return p.x < center.x + radius &&
               p.x > center.x - radius &&
               p.y < center.y + radius &&
               p.y > center.y - radius;
    }
};

// each QuadTree contains one point, initializes to square of radius 2 centered at origin
struct QuadTree {
    // true if node is a leaf, aka no children
    bool is_leaf;
    // four children
    QuadTree* nW;
    QuadTree* nE;
    QuadTree* sW;
    QuadTree* sE;

    // bounding box to represent the boundaries of this quad tree
    Box boundary;

    // true if there are no points within this tree, aka no centroid
    bool is_empty;
    // centroid, contains coordinate and phase
    Point centroid;

    // "mass", aka how many points
    int mass;

    // constructor
    QuadTree(Box b = Box()):
        is_leaf(true), is_empty(true), nW(NULL), nE(NULL), sW(NULL), sE(NULL),
        centroid(Point()), mass(0), boundary(b) {}

    // create four children
    void subdivide(){
        double subradius = boundary.radius / 2;

        Point qcenter = Point(boundary.center.x - subradius, boundary.center.y + subradius);
        nW = new QuadTree(Box(qcenter, subradius));

        qcenter = Point(boundary.center.x + subradius, boundary.center.y + subradius);
        nE = new QuadTree(Box(qcenter, subradius));

        qcenter = Point(boundary.center.x - subradius, boundary.center.y - subradius);
        sW = new QuadTree(Box(qcenter, subradius));

        qcenter = Point(boundary.center.x + subradius, boundary.center.y - subradius);
        sE = new QuadTree(Box(qcenter, subradius));

        is_leaf = false;
    }

    // insert point into quadtree, updating its centroid in the process
    void insert(Point p) {
        // ignore objects that are not in current bounds, this should never happen
        if (!boundary.contains(p)) throw;

        // if there is space in this quad tree, add point
        if (is_empty) {
            centroid = p;
            mass++;
            is_empty = false;
            return;
        }

        // this tree already has a point, aka already has a centroid

        // if no children, subdivide
        if (is_leaf) {
            subdivide();
            if (nW->boundary.contains(centroid)) nW->insert(centroid);
            else if (nE->boundary.contains(centroid)) nE->insert(centroid);
            else if (sE->boundary.contains(centroid)) sE->insert(centroid);
            else if (sW->boundary.contains(centroid)) sW->insert(centroid);
            else throw;
        }

        // update current tree centroid
        // weighted averages of coordinate and phase
        int m_new = mass + 1;
        double x_bar = (mass * centroid.x + p.x) / m_new;
        double y_bar = (mass * centroid.y + p.y) / m_new;
        double phase_bar = atan2((mass * sin(centroid.phase) + sin(p.phase)) / m_new,
                                 (mass * cos(centroid.phase) + cos(p.phase)) / m_new);
        centroid = Point(x_bar, y_bar, phase_bar);
        mass = m_new;

        // find new children that will eventually accept this point
        if (nW->boundary.contains(p)) nW->insert(p);
        else if (nE->boundary.contains(p)) nE->insert(p);
        else if (sE->boundary.contains(p)) sE->insert(p);
        else if (sW->boundary.contains(p)) sW->insert(p);
    }

    std::vector<double> get_centroids(double x, double y, double theta) {
        std::vector<double> out;

        if (is_empty){
            return out;
        }

        double cx = boundary.center.x;
        double cy = boundary.center.y;
        double cw = 2 * boundary.radius;

        if (is_leaf ||
            cw / (sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy))) < theta) {

            if (centroid.x == x && centroid.y == y) return out;

            out.push_back(centroid.x);
            out.push_back(centroid.y);
            out.push_back(centroid.phase);
            out.push_back(mass);

            return out;
        }

        std::vector<double> v1, v2, v3, v4;
        v1 = nW->get_centroids(x, y, theta);
        v2 = nE->get_centroids(x, y, theta);
        v3 = sE->get_centroids(x, y, theta);
        v4 = sW->get_centroids(x, y, theta);

        out.insert(out.end(), v1.begin(), v1.end());
        out.insert(out.end(), v2.begin(), v2.end());
        out.insert(out.end(), v3.begin(), v3.end());
        out.insert(out.end(), v4.begin(), v4.end());

        return out;
    }

    // destructor
    ~QuadTree() {
        delete nW; delete nE; delete sW; delete sE;
    }
};
