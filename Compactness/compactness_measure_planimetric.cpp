#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

//define CGAL_CLASSIFICATION_VERBOSE true
/*
#include <boost/timer/progress_display.hpp>
#ifdef _USE_OPENMP
#include <omp.h>
std::cout << "omp.h"
#endif
*/
#include <omp.h>
#include <cstdlib>
#include <memory> // std::make_unique
#include <fstream>
#include <iostream>
#include <string>
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
using Planimetric_grid = Classification::Planimetric_grid<Kernel, Point_set, Pmap>;

/*
#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
namespace AW3 = CGAL::Alpha_wraps_3;
using Mesh = CGAL::Surface_mesh<Point>; */

#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
namespace AW3 = CGAL::Alpha_wraps_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_Alpha;
typedef Kernel_Alpha::Point_3 Point_Alpha;
using Point_container = std::vector<Point_Alpha>;
using Mesh = CGAL::Surface_mesh<Point_Alpha>;

typedef boost::graph_traits<Mesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor            edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

// User-defined feature
#include "compactness.h"

int main (int argc, char** argv)
{
  if (argc == 1) {
      std::cout << "Usage: " << argv[0] << " filename scale n_scales alpha offset" << std::endl;
  }
  const std::string filename = (argc > 1) ? argv[1] : "C:/Users/hmalet/Documents/Myproject/Chambon_large_displacements/Chambon_Scan_Riegl_20220630 - small.ply";
  std::ifstream in (filename.c_str(), std::ios::binary);
  Point_set pts;
  std::cerr << "Reading input" << std::endl;
  in >> pts;
  if (pts.number_of_points() == 0) {
      std::cerr << "Error: no vertices found." << std::endl;
      return EXIT_FAILURE;
  }

std::cout << "No. of vertices: " << pts.number_of_points() << std::endl;
  
  std::cout << "Properties found:" << std::endl;
  for (auto prop : pts.properties_and_types()) {
      std::cout << "  " << prop.first << std::endl;
  }

  // radius_size = 0.6; relative_alpha = 0.1; relative_offset = 2.0;

  const float radius_size = (argc > 2) ? atof(argv[2]) : 0.3f; // specify radius of neighborhoods (default: 60cm, MAC suggests: 10cm)
  //const float voxel_size = radius_size / 3.f; // re-scale for CGAL's feature generator
  const float voxel_size = (radius_size*5) / 3.f; // re-scale for CGAL's feature generator (multiply by 5 to not have problems with too fine of a grid)
 
  // wrap surface
  // Compute the alpha and offset values
  const double relative_alpha = (argc > 4) ? std::stod(argv[4]) : 0.15f;//0.2;// 10. 
  const double relative_offset = (argc > 5) ? std::stod(argv[5]) : 2.f;//2.f;// 300.;
  std::cout << "relative alpha = " << relative_alpha << " relative offset = " << relative_offset << std::endl;
  double alpha = radius_size / relative_alpha; // bbox / relative_alpha;
  double offset = radius_size / relative_offset; // bbox / relative_offset
  std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;

  CGAL::Real_timer t;
  t.start();

  // convert to a kernel that is more stable for Alpha Wrap
  Point_container points;
  for (auto& point : pts.points()) {
      Point_Alpha pt(point.x(), point.y(), point.z());
      points.push_back(pt);
  }

  // construct the wrap
  Mesh wrap;
  CGAL::alpha_wrap_3(points, alpha, offset, wrap);
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, " << std::endl;
  
  Mesh::Property_map<edge_descriptor, bool> is_constrained_map = wrap.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;

  CGAL::IO::write_polygon_mesh("wrap.ply", wrap, CGAL::parameters::stream_precision(17));
  std::cout << "Wrap saved" << std::endl;

  t.stop();
  std::cout << "Took " << t.time() << " s" << std::endl;


  Feature_set features;
  std::cerr << "Initialising feature generator...";

  //const float radius_size = (argc > 2) ? atof(argv[2]) : 0.1f; // specify radius of neighborhoods (default: 10cm)
  t.reset(); t.start();
  int n_scales = (argc > 3) ? std::stod(argv[3]) : 2; // normally 5
  Feature_generator generator(pts, pts.point_map(), n_scales, voxel_size);
  t.stop();
  std::cout << "done in " << t.time() << " second(s)" << std::endl;

  std::cout << "Neighbourhood radii: " << generator.radius_neighbors(0);
  for (std::size_t i = 1; i < generator.number_of_scales(); ++i)
      std::cout << ", " << generator.radius_neighbors(i);
  std::cout << std::endl;
  
  std::cout << "Computing our features..." << std::endl;
  t.reset();
  t.start();
  auto bbox = CGAL::bounding_box(CGAL::make_transform_iterator_from_property_map(pts.begin(), pts.point_map()),
                                 CGAL::make_transform_iterator_from_property_map(pts.end(), pts.point_map()));
  using MyFeature = CGAL::Classification::Feature::My_feature<Kernel, Point_set, Pmap>;
 //
  for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
      Planimetric_grid grid(pts, pts.point_map(), bbox, generator.grid_resolution(i));
      features.add_with_scale_id<MyFeature>(i, pts, pts.point_map(), grid, generator.radius_neighbors(i), wrap);
  }
  t.stop();
  std::cout << "  done in " << t.time() << " second(s)" << std::endl;

  // output results
  for (auto& feature : features) {
      Fmap fmap = pts.add_property_map<float>("scalar_" + feature->name(), 0).first;
      for (std::size_t i = 0; i < pts.size(); ++i) {
          if (i < 5) std::cout << feature->name() << "[" << i << "] " << feature->value(i) << std::endl;
          fmap[i] = feature->value(i);
      }
  }

  std::cout << "Writing output..." << std::endl;
  std::ofstream f("output.ply");
  f.precision(18);
  f << pts;
  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}