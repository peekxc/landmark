////////////////////////////////////////////////////////////////////////////////
// Maxmin-based Landmark Procedure
// Authors: Matt Piekenbrock, Jason Cory Brunson, Yara Skaf
// Description: Calculate a landmark set using the maxmin procedure.
////////////////////////////////////////////////////////////////////////////////

#include "rankdistance.h"
using namespace Rcpp;
using namespace std;
using std::size_t;

// define a tolerance for distances to avoid floating point rounding errors
#define EPS_TOL 0.001

inline double max_dist(const NumericVector& x, const NumericVector& y){
  return(max(abs(x-y)));
}

inline double man_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = abs(x-y);
  return(sum(diff));
}

inline double sq_dist(const NumericVector& x, const NumericVector& y){
  NumericVector diff = x-y;
  return(sqrt(sum(pow(diff, 2.0))));
}

// Description:
//   Maxmin procedure to choose 'num' landmarks for balls of fixed radius
//   (supports euclidean distance only)
//
// Parameters:
//   - x : a data matrix
//   - num : desired number of landmark points, or number of sets, in a
//     ball cover (should be a positive integer)
//   - radius : desired radius of a cover set (should be a positive real number)
//   - seed_index : index of the first landmark used to seed the algorithm
//   - cover : boolean specifying whether to return cover sets in addition to
//     the landmark points
//
// [[Rcpp::export]]
List landmarks_maxmin_cpp(const NumericMatrix& x, int num = 0, float radius = -1, const int seed_index = 1, const bool cover = false) {
    int num_pts = x.nrow();

    // error handling
    if(radius < 0 && radius != -1){stop("Parameter 'radius' must be a positive number.");}
    if(num < 0){stop("Parameter 'num' must be >= 1.");}
    if(seed_index < 1 || seed_index > num_pts){stop("Parameter 'seed_index' must be >=1 and <= number of data points.");}

    // additional parameter handling
    if(num > num_pts){
        warning("Warning: parameter 'num' was > max allowable value. Setting num = number of data points.");
        num = num_pts;
    }
    if(num == 0 && radius == -1){num = std::min(num_pts,24);} // no parameters passed -> default behavior
    if(radius == -1){radius = FLT_MAX;} // Rcpp does not allow use of c++ constants (e.g. FLT_MAX) in parameters

    // store indices and values of X\L
    map<int, vector<double>> pts_left;
    for(int i = 0; i < num_pts; i++){
        NumericVector vec = x.row(i);
        vector<double> point(vec.begin(),vec.end());
        pts_left.emplace(i, point);
    }

    // remove seed landmark (and any duplicates) from X\L
    vector<double> seed_val = pts_left.at(seed_index-1);
    for (auto pt = pts_left.cbegin(), next_pt = pt; pt != pts_left.cend(); pt = next_pt){
      ++next_pt;
      if (pt->second == seed_val){ pts_left.erase(pt); }
    }

    // store indices and values of landmark set L
    map<int, vector<double>> landmarks;
    vector<int> ordered_landmarks; // keep track of order landmarks are added
    landmarks.emplace(seed_index-1, seed_val);
    ordered_landmarks.push_back(seed_index-1);

    // compute remaining landmarks
    while(true){
        map<int, vector<double>> maxmin;
        double d_max = 0;

        // find max(d(x,L)) for x in X
        for(const auto& pt : pts_left){
            double d_min = DBL_MAX;

            // find min(d(x,l)) for l in L
            for(const auto& l : landmarks){
                double d = dist_euc(l.second, pt.second);
                if(d < d_min){ d_min = d; }
            }

            // d_min is equal to the old max -> add this point to maxmin
            if(d_min == d_max){ maxmin.insert(pt); }

            // we have a new max -> clear out maxmin and add this point instead
            if(d_min > d_max){
                d_max = d_min;
                maxmin.clear();
                maxmin.insert(pt);
            }
        }
        // done if farthest point is within radius and we have enough landmarks
        if(d_max <= radius && landmarks.size() >= num){
            if(radius == FLT_MAX){ radius = d_max; }
            break;
        }

        // otherwise add new max to L and remove from X\L
        pair<int, vector<double>> l_i = make_pair(maxmin.begin()->first, maxmin.begin()->second);
        ordered_landmarks.push_back(l_i.first);
        landmarks.insert(l_i);

        // remove all points at this center & exit if all points are covered
        for(const auto& pt : maxmin){ if(pt.second == l_i.second){ pts_left.erase(pt.first); } }
        if(pts_left.size() <= 0){
            // all points coincide with a point in L, so ball radius should be 0
            if(radius == FLT_MAX){ radius = 0; }
            break;
        }
    }

    List ret;
    IntegerVector landmarks_R = wrap(ordered_landmarks); // wrap into R data type
    ret.push_back(landmarks_R+1);

    if(cover == true){
        List cover_sets;
        for(const auto& l : ordered_landmarks){
            vector<double> l_val = landmarks.at(l);
            IntegerVector ball;
            for(int i = 0; i < num_pts; i++){
                // if pt is within the radius (+ tolerance), add it to the ball
                NumericVector vec = x.row(i);
                vector<double> pt(vec.begin(),vec.end());
                if(dist_euc(l_val, pt) <= radius+EPS_TOL){ ball.push_back(i); }
            }
            cover_sets.push_back(ball+1);
        }
        ret.push_back(cover_sets);
    }
    return(ret);
}

// a Distance Function here just returns the distance between two points, given their *indices*
using DistFunction = typename std::function<double(size_t, size_t)>;

// dist_f  := distance function between two (indexed) points
// n_pts   := number of points in the data set
// eps     := distance threshold used as a stopping criterion for the maxmin procedure
// n       := cardinality threshold used as a stopping criterion for the maxmin procedure
// metric  := metric to use. If 0, uses `dist_f`, otherwise picks one of the available metrics.
// seed    := initial point (default is point at index 0)
// pick    := criterion to break ties. Possible values include 0 (first), 1 (random), or 2 (last).
// cover   := whether to report set membership for each point
List maxmin_f(DistFunction dist_f, const size_t n_pts,
              const double eps, const size_t n,
              const size_t seed = 0, const size_t pick = 0,
              const bool cover = false) {
  if (eps == -1.0 && n == 0){ stop("Must supply either positive 'eps' or positive 'n'."); }
  if (pick > 2){ stop("tiebreaker 'pick' choice must be in { 0, 1, 2 }."); }
  if (seed >= n_pts){ stop("Invalid seed index given."); }

  // Make a function that acts as a sentinel
  enum CRITERION { NUM, EPS, NUM_OR_EPS }; // These are the only ones that make sense
  const CRITERION stopping_criterion = (eps == -1.0) ? NUM : ((n == 0) ? EPS : NUM_OR_EPS);
  const auto is_finished = [stopping_criterion, eps, n](size_t n_landmarks, double c_eps){
    switch(stopping_criterion){
      case NUM: return(n_landmarks >= n);
      case EPS: return(c_eps <= eps);
      case NUM_OR_EPS: return(n_landmarks >= n || c_eps <= eps);
    }
  };

  // Indices of possible candidate landmarks
  vector< size_t > candidate_pts(n_pts);
  std::iota(begin(candidate_pts), end(candidate_pts), 0);

  // Choose the initial landmark
  vector< size_t > lm = vector< size_t >();
  if (n != 0){ lm.reserve(n); }
  lm.push_back(seed);

  // Remove initial seed from candidates
  candidate_pts.erase(begin(candidate_pts) + seed);

  // Indices assign each point to its closest landmark index
  vector< size_t > closest_landmark(n_pts, seed);

  // Preallocate distance vector for landmarks; one for each point
  std::vector< double > lm_dist(n_pts, std::numeric_limits<double>::infinity());

  // Generate the landmarks
  bool stop_reached = false;
  while (!stop_reached){
    // Inductively, replace point-to-landmark distances if lower than previously computed
    const size_t c_lm = lm.back();

    // Update non-landmark points with distance to nearest landmark
    for (auto idx: candidate_pts){
      double c_dist = dist_f(c_lm, idx);
      if (c_dist < lm_dist[idx]){
        lm_dist[idx] = c_dist;                  // update minimum landmark distance
        closest_landmark[idx] = lm.size() - 1;  // mark which landmark was closest, by index
      }
    }

    // Of the remaining candidate points, find the one with the maximum landmark distance
    auto max_landmark = std::max_element(begin(candidate_pts), end(candidate_pts), [&lm_dist](size_t ii, size_t jj){
      return lm_dist[ii] < lm_dist[jj];
    });

    // If not greedily picking the first candidate point, partition the candidate points, then use corresponding strategy
    if (pick > 0 && max_landmark != end(candidate_pts)){
      double max_lm_dist = lm_dist[(*max_landmark)];
      auto it = std::partition(begin(candidate_pts), end(candidate_pts), [max_lm_dist, &lm_dist](size_t j){
        return lm_dist[j] == max_lm_dist;
      });
      // Pick == 1 => the last candidate
      if (pick == 1){
        max_landmark = it != begin(candidate_pts) ? std::prev(it) : begin(candidate_pts);
      } else {
      // Pick > 1 => random candidate with max distance
        const size_t n_cand = std::distance(begin(candidate_pts), it);
        max_landmark = (n_cand == 1) ? begin(candidate_pts) : begin(candidate_pts) + (rand() % n_cand);
      }
    }

    // If the iterator is valid, we have a new landmark, otherwise we're finished
    if (max_landmark != end(candidate_pts)){
      // Rprintf("max lm dist: %g\n", lm_dist[(*max_landmark)]);
      lm.push_back(*max_landmark);
      stop_reached = is_finished(lm.size(), lm_dist[(*max_landmark)]);
      candidate_pts.erase(max_landmark);
    } else {
      stop_reached = true;
    }
  } // while(!finished())

  // Prepare return
  IntegerVector lm_int = wrap(lm);
  List ret = List::create(_["landmarks"] = lm_int+1, _["cover"] = R_NilValue);

  // If requested, collect points into a cover
  if (cover){
    vector< IntegerVector > cover_sets(lm.size());
    for (size_t i = 0; i < n_pts; ++i){
      cover_sets.at(closest_landmark[i]).push_back(i+1); // +1 for 1-based indices
    }
    ret["cover"] = wrap(cover_sets);
    ret.attr("radius") = lm_dist[lm.back()]; // capture radius needed to make the cover
  }
  return(ret);
}

// Point cloud wrapper - See 'maxmin_f' below for implementation
// [[Rcpp::export]]
List maxmin_pc(const NumericMatrix& x, const double eps, const size_t n, Function dist_f,
               const size_t metric = 1, const size_t seed = 0, const size_t pick = 0,
               const bool cover = false){
  const size_t n_pts = x.nrow();
  if (seed < 0 || seed >= n_pts){ stop("Invalid seed point."); }

  // Choose the distance function
  DistFunction dist;
  if (metric == 0){
    dist = [&x, dist_f](size_t i, size_t j) { return as< double >(dist_f(x.row(i), x.row(j))); };
  } else if (metric == 1){
    dist = [&x](size_t i, size_t j) { return sq_dist(x.row(i), x.row(j)); };
  } else if (metric == 2){
    dist = [&x](size_t i, size_t j) { return man_dist(x.row(i), x.row(j)); };
  } else if (metric == 3){
    dist = [&x](size_t i, size_t j) { return max_dist(x.row(i), x.row(j)); };
  } else {
    stop("Invalid distance metric option chosen.");
  }

  // Call the generalized procedure
  return maxmin_f(dist, n_pts, eps, n, seed, pick, cover);
}

// Converts (i,j) indices in the range { 0, 1, ..., n - 1 } to its 0-based position
// in a lexicographical ordering of the (n choose 2) combinations.
constexpr size_t to_nat_2(size_t i, size_t j, size_t n) noexcept {
  return i < j ? (n*i - i*(i+1)/2 + j - i - 1) : (n*j - j*(j+1)/2 + i - j - 1);
}

// 'dist' object wrapper - See 'maxmin_f' below for implementation
// [[Rcpp::export]]
List maxmin_dist(const NumericVector& x, const size_t n_pts,
               const double eps, const size_t n,
               const size_t seed = 0, const size_t pick = 0,
               const bool cover = false){
  if (seed < 0 || seed >= n_pts){ stop("Invalid seed point."); }

  // Parameterize the distance function
  DistFunction dist = [&x, n_pts](size_t i, size_t j) -> double {
    return x[to_nat_2(i,j,n_pts)];
  };

  // Call the generalized procedure
  return maxmin_f(dist, n_pts, eps, n, seed, pick, cover);
}



