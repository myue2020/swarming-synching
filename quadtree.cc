#include <cmath>
#include <vector>
#include "quadtree.hh"

double LEAF = 1000000.;

// constructor
void QuadTree::QuadTree(Box b = Box())
{
    nW = NULL;
    nE = NULL;
    sW = NULL;
    sE = NULL;
    centroid = Point(LEAF, 0., 0.);
    mass = 0;
    boundary = b;
}

// destructor
void QuadTree::~QuadTree()
{
    delete nW; delete nE; delete sW; delete sE;
}

// create four children
void QuadTree::subdivide()
{
    double subradius = boundary.radius / 2;

    Point qcenter = Point(boundary.center.x - subradius, boundary.center.y + subradius);
    nW = new QuadTree(Box(qcenter, subradius));

    qcenter = Point(boundary.center.x + subradius, boundary.center.y + subradius);
    nE = new QuadTree(Box(qcenter, subradius));

    qcenter = Point(boundary.center.x - subradius, boundary.center.y - subradius);
    sW = new QuadTree(Box(qcenter, subradius));

    qcenter = Point(boundary.center.x + subradius, boundary.center.y - subradius);
    sE = new QuadTree(Box(qcenter, subradius));
}

// insert a point into the QuadTree, subdividing if needed
bool QuadTree::insert(Point p)
{
    // ignore objects that are not in current bounds, this should never happen
    if (!boundary.contains(p)) return false;

    // if there is space in this quad tree and no subdivisions
    if (centroid.x == LEAF){
        centroid = p;
        mass++;
        return true;
    }

    // subdivide
    if (nW == NULL){
        subdivide();
        if (nW->boundary.contains(centroid)) nW->insert(centroid);
        if (nE->boundary.contains(centroid)) nE->insert(centroid);
        if (sE->boundary.contains(centroid)) sE->insert(centroid);
        if (sW->boundary.contains(centroid)) sW->insert(centroid);
    }

    // find new children that will eventually accept this point
    if (nW->insert(p)) {mass++; return true;}
    if (nE->insert(p)) {mass++; return true;}
    if (sW->insert(p)) {mass++; return true;}
    if (sE->insert(p)) {mass++; return true;}

    return false;
}

// update centroid's coordinates and phase by traversing tree
void QuadTree::update_centroid()
{
    // only traverse further if there are children
    if (nW != NULL){
        Point p1, p2, p3, p4, avg;
        int m1=0, m2=0, m3=0, m4=0, m;
        double x_bar, y_bar, phase_bar;

        // if child has centroid, then update and retrieve it
        if (nW->centroid.x == LEAF){
            nW->update_centroid();
            p1 = nW->centroid;
            m1 = nW->mass;
        }
        if (nE->centroid.x == LEAF){
            nE->update_centroid();
            p2 = nE->centroid;
            m2 = nE->mass;
        }
        if (sE->centroid.x == LEAF){
            sE->update_centroid();
            p3 = sE->centroid;
            m3 = sE->mass;
        }
        if (sW->centroid.x == LEAF){
            sW->update_centroid();
            p4 = sW->centroid;
            m4 = sW->mass;
        }

        // weighted averages of coordinate and phase
        m = m1 + m2 + m3 + m4;
        x_bar = (m1*p1.x+m2*p2.x+m3*p3.x+m4*p4.x)/m;
        y_bar = (m1*p1.y+m2*p2.y+m3*p3.y+m4*p4.y)/m;
        phase_bar = atan2((m1*sin(p1.phase)+m2*sin(p2.phase)+m3*sin(p3.phase)+m4*sin(p4.phase))/m,
                          (m1*cos(p1.phase)+m2*cos(p2.phase)+m3*cos(p3.phase)+m4*cos(p4.phase))/m);
        avg = Point(x_bar, y_bar, phase_bar);
        centroid = avg;
    }
}

std::vector<double> QuadTree::get_centroids(double x, double y, double theta)
{
    std::vector<double> out;

    if (centroid.x == LEAF){
        return out;
    }

    double cx = boundary.center.x;
    double cy = boundary.center.y;
    double cw = 2 * boundary.radius;

    if (nW == NULL || cw/(sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy)))<theta){
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
