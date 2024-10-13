#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
//#include <boost/timer/progress_display.hpp>

// Classification
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Real_timer.h>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef Point_set::Point_map Pmap;
typedef Point_set::Property_map<int> Imap;
typedef Point_set::Property_map<float> Fmap;
typedef Point_set::Property_map<unsigned char> UCmap;
namespace Classification = CGAL::Classification;
typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;
typedef Classification::Point_set_feature_generator<Kernel, Point_set, Pmap>    Feature_generator;

using Neighborhood = Classification::Point_set_neighborhood<Kernel, Point_set, Pmap>;
using Neighbor_query = typename Neighborhood::Sphere_neighbor_query;

// Alpha wrap
#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
namespace AW3 = CGAL::Alpha_wraps_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_Alpha;
typedef Kernel_Alpha::Point_3 Point_Alpha;
using Point_container = std::vector<Point_Alpha>;
using Mesh = CGAL::Surface_mesh<Point_Alpha>;

// Sphere generation
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Aff_transformation_3.h>
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT(*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
FT sphere_function(Point_3 p) {
    const FT x2 = p.x() * p.x(), y2 = p.y() * p.y(), z2 = p.z() * p.z();
    return x2 + y2 + z2 - 1;
}

// Mesh intersection
#include <CGAL/Polygon_mesh_processing/corefinement.h>
typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor            vertex_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;


    // User-defined feature
    template <typename GeomTraits, typename PointRange, typename PointMap>
    class My_feature : public CGAL::Classification::Feature_base

    {

        using FloatMap = typename PointRange::template Property_map<float>;
        PointMap point_map;
        std::vector<typename FloatMap::value_type> values;

    public:
        My_feature(const PointRange& input, PointMap point_map, const float feature_scale, const Mesh& wrap) : point_map(point_map)//, neighbor_query(neighbor_query)
        {
            //const double feature_scale = (argc > 4) ? std::stod(argv[4]) : alpha * 10.0f; // feature scale (radius of sphere)
            this->set_name("my_feature");

            values.resize(input.size());

            Mesh::Property_map<edge_descriptor, bool> is_constrained_map = wrap.property_map<edge_descriptor, bool>("e:is_constrained").first;

            Mesh sm;
            // Create Sphere
            Tr tr;            // 3D-Delaunay triangulation
            C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation
            // defining the surface
            Surface_3 surface(sphere_function,             // pointer to function
                Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
    // defining meshing criteria
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                0.1,  // radius bound
                0.1); // distance bound
    // meshing surface
            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

            CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

            std::cout << feature_scale << std::endl;

            // carve out part of wrap around each point
            for (std::size_t i = 0; i < input.size(); i++) {
                if (i % int(input.size() / 100.f) == 0) {
                    std::cout << float(i) / input.size() * 100.f << "%" << std::endl;
                }
                // ***************************
                // RMD To speed up the code (for testing purposes) you might want to only compute this descriptor on a subset of the points initially.
                //if (i > 100) break;
                // ***************************
                const Point& point = input.point(i);
                Mesh smi = sm; // copy sphere mesh

                // Translate and scale sphere mesh
                for (Mesh::Vertex_index p : smi.vertices()) {
                    Point_3& a = smi.point(p);
                    smi.point(p) = Point_3(point.x() + a.x() * feature_scale, point.y() + a.y() * feature_scale, point.z() + a.z() * feature_scale);
                }

                // Intersect wrap and sphere
                Mesh intersected_mesh = wrap;
                bool valid_difference = PMP::corefine_and_compute_intersection(intersected_mesh,
                    smi,
                    intersected_mesh,
                    params::default_values(),
                    params::default_values(),
                    params::edge_is_constrained_map(intersected_mesh.property_map<edge_descriptor, bool>("e:is_constrained").first));

                if (i == 0)
                    CGAL::IO::write_polygon_mesh("intersected_mesh.ply", intersected_mesh, CGAL::parameters::stream_precision(17));


                // ***************************
                // RMD you want to replace this code for computing "Betti numbers" on the intersected mesh with code to compute the surface area and volume of intersected_mesh.
                // ***************************


                //compute volume 

                double vol = 0;

                for (Mesh::Face_range::iterator it = intersected_mesh.faces_begin(); it != intersected_mesh.faces_end(); ++it) {
                    halfedge_descriptor hd = intersected_mesh.halfedge(*it);

                    Point_container P;
                    do {
                        const Point_3& p = intersected_mesh.point(intersected_mesh.source(hd));
                        P.push_back(p);
                        hd = intersected_mesh.next(hd);
                    } while (hd != intersected_mesh.halfedge(*it));

                    // Determinant
                    double v321 = P.at(2).x() * P.at(1).y() * P.at(0).z();
                    double v231 = P.at(1).x() * P.at(2).y() * P.at(0).z();
                    double v312 = P.at(2).x() * P.at(0).y() * P.at(1).z();
                    double v132 = P.at(0).x() * P.at(2).y() * P.at(1).z();
                    double v213 = P.at(1).x() * P.at(0).y() * P.at(2).z();
                    double v123 = P.at(0).x() * P.at(1).y() * P.at(2).z();

                    double det = -v321 + v231 + v312 - v132 - v213 + v123;
                    vol += det;

                }
                vol /= (double)6.0f;
                if (i == 0)
                    std::cout << "Volume :" << vol << std::endl;

                //compute surface area 

                double surface = 0;

                for (Mesh::Face_range::iterator it = intersected_mesh.faces_begin(); it != intersected_mesh.faces_end(); ++it) {
                    halfedge_descriptor hd = intersected_mesh.halfedge(*it);

                    Point_container P;
                    do {
                        const Point_3& p = intersected_mesh.point(intersected_mesh.source(hd));
                        P.push_back(p);
                        hd = intersected_mesh.next(hd);
                    } while (hd != intersected_mesh.halfedge(*it));
                    
                    double area = 0.5 * sqrt(std::pow((P.at(1).y() - P.at(0).y()) * (P.at(2).z() - P.at(0).z()) - (P.at(1).z() - P.at(0).z()) * (P.at(2).y() - P.at(0).y()), 2) +
                        (std::pow((P.at(1).z() - P.at(0).z()) * (P.at(2).x() - P.at(0).x()) - (P.at(1).x() - P.at(0).x()) * (P.at(2).z() - P.at(0).z()), 2)
                            + (std::pow((P.at(1).x() - P.at(0).x()) * (P.at(2).y() - P.at(0).y()) - (P.at(1).y() - P.at(0).y()) * (P.at(2).x() - P.at(0).x()), 2))));

                    surface += std::abs(area);
                }
                if (i == 0)
                    std::cout << "Surface area :" << surface << std::endl;

                double compactness = (3 * vol) / (feature_scale * surface); 
                values[i] = compactness;
                if (i == 0)
                    std::cout << "compactness :" << compactness << std::endl;
            }
        }
        float value(std::size_t pt_index) {
            return values[pt_index];
        }
    };

    
    int main(int argc, char** argv) {
        if (argc == 1) std::cout << "Usage: " << argv[0] << " filename scale alpha offset" << std::endl;

        // ***************************
        // RMD change the file name
        // ***************************
        std::string filename = (argc > 1) ? argv[1] : "C:/Users/hmalet/Documents/Myproject/Chambon_large_displacements/Chambon_Scan_Riegl_20220630 - small.ply";

        std::ifstream in(filename.c_str(), std::ios::binary);
        Point_set pts;
        std::cerr << "Reading input" << std::endl;
        in >> pts;
        if (pts.number_of_points() == 0) {
            std::cerr << "Error: no vertices found." << std::endl;
            return EXIT_FAILURE;

        }

        std::cout << "No. of vertices: " << pts.number_of_points() << std::endl;

        const float radius_size = (argc > 2) ? atof(argv[2]) : 0.1f; // specify radius of neighborhoods (default: 10cm)
        const float voxel_size = radius_size / 3.f; // re-scale for CGAL's feature generator

        // wrap surface
        // Compute the alpha and offset values
        const double relative_alpha = (argc > 2) ? std::stod(argv[3]) : 0.05;// 10.;
        const double relative_offset = (argc > 3) ? std::stod(argv[4]) : 2.f;// 300.;
        std::cout << "relative alpha = " << relative_alpha << " relative offset = " << relative_offset << std::endl;
        double alpha = radius_size / relative_alpha; // bbox / relative_alpha;
        double offset = radius_size / relative_offset; // bbox / relative_offset
        std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;
        alpha = 2.5; // 5.
        offset = 0.5;
        CGAL::Real_timer t;
        t.start();

        // convert to a kernel that is more stable for Alpha Wrap
        Point_container points;
        for (auto& point : pts.points()) {
            Point_Alpha pt(point.x(), point.y(), point.z());
            points.push_back(pt);
        }

        // construct the wrap
        std::cout << "Wrapping point cloud...";
        Mesh wrap;
        CGAL::alpha_wrap_3(points, alpha, offset, wrap);
        t.stop();
        std::cout << "done (" << t.time() << " s)" << std::endl;

        std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, " << std::endl;

        Mesh::Property_map<edge_descriptor, bool> is_constrained_map = wrap.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;

        CGAL::IO::write_polygon_mesh("wrap.ply", wrap, CGAL::parameters::stream_precision(17));
        std::cout << "Wrap saved" << std::endl;



        std::cout << "Computing our features...";
        t.reset();
        t.start();
        std::unique_ptr<Neighborhood> neighborhood = std::make_unique<Neighborhood>(pts, pts.point_map());
        using MyFeature = My_feature<Kernel, Point_set, Pmap>;
        MyFeature my_feature(pts, pts.point_map(), radius_size, wrap);
        t.stop();
        std::cout << "done (" << t.time() << " s)" << std::endl;

        // output results
        Fmap compactness = pts.add_property_map<float>("scalar_Compactness", 0).first;
        for (std::size_t i = 0; i < pts.size(); ++i)
        {
            if (i < 5) std::cout << "compactness[" << i << "] " << my_feature.value(i) << std::endl;
            compactness[i] = my_feature.value(i);
        }

        std::cout << "Writing output..." << std::endl;
        std::ofstream f("output.ply");
        f.precision(18);
        f << pts;
        std::cerr << "All done" << std::endl;

        return EXIT_SUCCESS;

    }
