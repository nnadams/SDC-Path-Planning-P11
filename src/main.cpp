#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

#define MAX_PATH_SIZE 50

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

enum State { KeepLane, PrepareLaneChange, Wait };
inline int lane2d(int lane) { return 2+4*lane; }
inline int d2lane(double d) { return (int)round((d-2)/4); }
inline bool isInLane(int lane, double car_d, double other_d) { return (other_d < (lane2d(lane)+2)) && (other_d > (lane2d(lane)-2)); }
bool isInAdjacentLane(int lane, double car_d, double other_d) {
  if (lane == 0 || lane == 2)
    return isInLane(1, car_d, other_d);
  else if (lane == 1)
    return isInLane(0, car_d, other_d) || isInLane(2, car_d, other_d);
  else
    return false;
}
bool isBlocking(int lane, double car_s, double car_d, double other_s, double other_d) { return (other_s < (car_s+30)) && (other_s > (car_s-8)); }

void process_state(State& state, int& car_lane, double& ref_vel, const double car_s, const double car_d, const int prev_size, const vector<vector<double>>& sensors) {
  static int refractory = 0;
  bool car_ahead_close = false;
  bool car_ahead_slow = false;
  float lane_cost[3] = {0,0,0};
  State next_state = state;

  cout << "State: " << state << endl;
  for (int i = 0; i < sensors.size(); i++) {
    double x = sensors[i][1];
    double y = sensors[i][2];
    double vx = sensors[i][3];
    double vy = sensors[i][4];
    double s = sensors[i][5];
    double d = sensors[i][6];
    int lane = d2lane(d);
    double speed = sqrt(vx*vx+vy*vy);

    // Attempt to handle latency
    s += ((double) prev_size*.02*speed);

    cout << "Car #" << sensors[i][0] << " in lane " << lane << " at (" << s << "," << d << ") going " << speed;

    // Set flags and cost for each lane
    if (isInLane(car_lane, car_d, d)) {
        if ((s > car_s) && ((s-car_s) < 30)) {
          car_ahead_close = true;
        }
        if (state == KeepLane) {
          if ((s > car_s) && ((s-car_s) < 50)) {
            next_state = PrepareLaneChange;
          }
        }
    }
    if (isBlocking(car_lane, car_s, car_d, s, d)) {
      lane_cost[lane] += 1000;
      cout << " +1000 Lane #" << lane << " is blocked!";
    }
    else if ((s-10) > car_s) {
      lane_cost[lane] += (75.0 / fabs(s - car_s)) + (3*speed);
      cout << " +" << lane_cost[lane];
    }
    cout << endl;
  }
  cout << endl << endl;

  // Slow down if close to the car ahead
  if (car_ahead_close) {
    ref_vel -= .25;
  }
  // Speed up if going too slow
  else if (ref_vel < 49.5) {
    ref_vel += .3;
  }

  cout << lane_cost[0] << ", " << lane_cost[1] << ", " << lane_cost[2] << endl;

  // If the cost of the current lane is high, try to move to an empty lane
  if (lane_cost[car_lane] > 35) {
    for (int l = 0; l < 3; l++) {
      if (lane_cost[l] == 0) next_state = PrepareLaneChange;
    }
  }

  if (state == Wait) next_state = KeepLane;
  else if (state == PrepareLaneChange && !refractory) {
    // We decided to try to change lanes, so pick a lane based on the cost
    int goal_lane = car_lane;
    int min_cars = INT_MAX;
    for (int l = 0; l < 3; l++) {
      if (l == car_lane) continue;
      if (lane_cost[l] < min_cars) {
        goal_lane = l;
        min_cars = lane_cost[l];
      }
    }

    // If the goal lane is 2 away just move to the center for now (less jerk)
    if (fabs(car_lane - goal_lane) > 1) goal_lane = 1;

    // If it's safe, set our new goal lane
    if (lane_cost[goal_lane] < 1000 && ref_vel > 25) {
      refractory = 50;
      car_lane = goal_lane;
      ref_vel -= .3;
    }
    next_state = KeepLane;
  }

  // Keep a refractory period so lane changes don't occur too often
  refractory = --refractory < 0 ? 0 : refractory;
  state = next_state;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  State car_state = KeepLane;
  int lane = 1; // vehicle starts in the middle lane
  double ref_vel = 0.0; // start velocity at 0 and slowly increase

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&ref_vel,&car_state](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            int prev_size = previous_path_x.size();

            if (prev_size > 0) {
              car_s = end_path_s;
            }

            cout << "Car at (" << car_s << "," << car_d << ") in #" << lane << endl;
            process_state(car_state, lane, ref_vel, car_s, car_s, prev_size, sensor_fusion);

            // The X,Y points on the road that the car will snap to
            vector<double> anchors_x;
            vector<double> anchors_y;

            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            if (prev_size < 2) { // If just starting only use the cars previous positions
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              anchors_x.push_back(prev_car_x);
              anchors_x.push_back(car_x);

              anchors_y.push_back(prev_car_y);
              anchors_y.push_back(car_y);
            }
            else { // Otherwise continue the path from the last iteration
              ref_x = previous_path_x[prev_size-1];
              ref_y = previous_path_y[prev_size-1];

              double ref_x_prev = previous_path_x[prev_size-2];
              double ref_y_prev = previous_path_y[prev_size-2];
              ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

              anchors_x.push_back(ref_x_prev);
              anchors_x.push_back(ref_x);

              anchors_y.push_back(ref_y_prev);
              anchors_y.push_back(ref_y);
            }

            // Use points further up the road
            vector<double> next_wp0 = getXY(car_s+30, lane2d(lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s+50, lane2d(lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s+90, lane2d(lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            anchors_x.push_back(next_wp0[0]);
            anchors_x.push_back(next_wp1[0]);
            anchors_x.push_back(next_wp2[0]);

            anchors_y.push_back(next_wp0[1]);
            anchors_y.push_back(next_wp1[1]);
            anchors_y.push_back(next_wp2[1]);

            // Transform anchor points to vehicle coordinates, with vehicle at origin
            for (int i = 0; i < anchors_x.size(); i++) {
              double shift_x = anchors_x[i] - ref_x;
              double shift_y = anchors_y[i] - ref_y;

              anchors_x[i] = (shift_x * cos(0-ref_yaw)) - (shift_y * sin(0-ref_yaw));
              anchors_y[i] = (shift_x * sin(0-ref_yaw)) + (shift_y * cos(0-ref_yaw));
            }

            // Use a spline to create a smooth path
            tk::spline s;
            s.set_points(anchors_x, anchors_y);

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            for (int i = 0; i < previous_path_x.size(); i++) {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            double target_x = 30.0;
            double target_y = s(target_x);
            double target_dist = sqrt((target_x*target_x)+(target_y*target_y));

            // Calculate each point (waypoints that are drawn in the simulator)
            double next_x = 0;
            int path_size = MAX_PATH_SIZE - previous_path_x.size();
            for (int i = 1; i <= path_size; i++) {
              double N = (target_dist / (.02*ref_vel/2.24));
              double x_point = next_x+target_x/N;
              double y_point = s(x_point);
              double x_ref = x_point;
              double y_ref = y_point;

              next_x = x_point;

              // Rotate and shift the points back to the global coordinates
              x_point = (x_ref * cos(ref_yaw)) - (y_ref * sin(ref_yaw));
              y_point = (x_ref * sin(ref_yaw)) + (y_ref * cos(ref_yaw));

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }

          	json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
