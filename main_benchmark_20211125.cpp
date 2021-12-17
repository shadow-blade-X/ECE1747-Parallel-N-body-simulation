//This is mainly based on the sequence given by https://conference.sdo.esoc.esa.int/proceedings/sdc7/paper/14/SDC7-paper14.pdf
//The timestep calculation is based on Gravitational N-Body Simulations, SVERRE J. AARSETH

#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <fstream>
#include <exception>
#include <sstream>
#include <functional>
#include <string>
#include "omp.h"


using namespace std;
std::hash<std::string> h;

const long double G=6.6743015e-11, PI=3.1415926535898;
vector <long double> add_init_position, add_init_velocity, body_values, output_vec,
					r_rel, v_rel, acting_values, a, a_dot, body_final,
					r_0,v_0,r_p,v_p,r_f,v_f,accel_output,a_0,a_0_dot,a_p,a_p_dot,
					a_i,a_i_dot,a_j,a_j_dot,r_i,v_i,r_j,v_j,a_ij,a_dot_ij,a_t_dot_ij,
					a_d_dot_ij,a_d_dot,a_t_dot,
					burn_vector, burn_vector_next, burn_ori, burn_ori_rate;
vector <string> load_results;
vector <string> body_string_ids;

long double add_mass, r_dot_v_relative, abs_r_rel, a_comp, a_dot_comp_1, a_dot_comp_2, abs_v_rel, comp_a, comp_b, comp_c, body_timestep, abs_a_t_dot, abs_a_d_dot, abs_a_dot, abs_a;
// string line, add_id, body_id, acting_id, burn_id, burn_body;
string line, add_id, burn_id, burn_body;
long double acting_id,body_id;

int itts, time_counter;
int burn_count = 0;
pair<string, vector<long double> > burn_values, burn_values_next;

class System {
	public:
		// map<string,vector<long double> > bodies, bodies_next;
        map<long double,vector<long double> > bodies, bodies_next;
		map <long double, string> body_id_double_to_string;
		// map<string, pair<string, vector<long double> > > burns;
		map<int, pair<long double, vector<long double> > > burns;
		long double accuracy, timestep, time, next_timestep;

		vector <vector<long double> > body_info;
		// vector <string> body_ids;
		vector <long double> body_ids;

		vector <vector<long double> > body_info_next;
		// vector <string> body_ids_next;
		vector <long double> body_ids_next;

		int output_rate;
		string info, debug_info;
        void SetAccuracy(int new_accuracy){
			accuracy=new_accuracy;
		}
		void SetStartTime(int start_time){
			time=start_time;
		}
		void SetOutToFile(const char * file_name, int rate){
			freopen(file_name,"w",stdout);
			output_rate = rate;
		}
		void CloseFile(){
			fclose(stdout);
		}

		void initialize(){
			// for (std::pair<std::string, vector <long double> > body_itterator : bodies){
			for (std::pair<long double, vector <long double> > body_itterator : bodies){
				body_ids.push_back(body_itterator.first);
				body_info.push_back(body_itterator.second);
			}
			body_info_next=body_info;

		}

		void AddBody (long double id, long double mass, vector<long double> init_position, vector<long double> init_velocity){
            bodies[id].clear();//so adding bodies with the same id twice doesn't break it
			bodies[id].push_back(mass);

			for( int i=0;i<init_position.size();i++){
				bodies[id].push_back(init_position.at(i));
			}

			for( int i=0;i<init_velocity.size();i++){
				bodies[id].push_back(init_velocity.at(i));
			}
		}

		// void AddBurn (string burn_id, string body_id, long double start_time, long double end_time, long double acceleration, vector<long double> orientation, vector<long double> orientation_rate){
		void AddBurn (int burn_id, long double body_id, long double start_time, long double end_time, long double acceleration, vector<long double> orientation, vector<long double> orientation_rate){
			//Gives burn_id:(body_id:[start,end,accel,orie_x,orie_y,orie_z,orie_ra_x...])
			burn_vector.clear();
			burn_vector.push_back(start_time);
			burn_vector.push_back(end_time);
			burn_vector.push_back(acceleration);
			for( int i=0;i<3;i++){
				burn_vector.push_back(orientation.at(i));
			}
			for( int i=0;i<3;i++){
				burn_vector.push_back(orientation_rate.at(i));
			}
			burns[burn_id] = make_pair(body_id,burn_vector);
		}
        void LoadFile(string filename){
			//Format is #body,id,mass,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z
			ifstream file(filename);
			if (file.is_open()) {
				string line;
				while (getline(file, line)) {
					if(line.compare(0,1,"#")==0){
						if(line.compare(1, 4,"Body")==0){
							load_results.clear();
							stringstream s_stream(line);
							while(s_stream.good()) {
							    string substr;
							    getline(s_stream, substr, ',');
							    load_results.push_back(substr);
							}
							add_init_position.clear();
							add_init_velocity.clear();

							add_id = load_results.at(1);

							std::cout << "This is add_id:" << add_id << "\n";
							std::cout << "This is add_id hash:" << h(add_id) << "\n";
							
							add_mass = stold(load_results.at(2));
							
							for( int i=0;i<3;i++){
								add_init_position.push_back(stold(load_results.at(3+i)));
								add_init_velocity.push_back(stold(load_results.at(6+i)));
							}
							body_id_double_to_string[h(add_id)]=add_id;
							body_string_ids.push_back(add_id);
							// AddBody(add_id, add_mass, add_init_position, add_init_velocity);
							AddBody(h(add_id), add_mass, add_init_position, add_init_velocity);

						}
						if(line.compare(1, 4,"Burn")==0){
							load_results.clear();
							stringstream s_stream(line);
							while(s_stream.good()) {
							    string substr;
							    getline(s_stream, substr, ',');
							    load_results.push_back(substr);
							}
							burn_ori.clear();burn_ori_rate.clear();
							for( int i=0;i<3;i++){
								burn_ori.push_back(stold(load_results.at(5+i)));
								burn_ori_rate.push_back(stold(load_results.at(8+i)));
							}
							burn_id = load_results.at(1);
							body_string_ids.push_back(burn_id);
							// AddBurn(to_string(burn_count),load_results.at(1),stold(load_results.at(2)),stold(load_results.at(3)),stold(load_results.at(4)),burn_ori,burn_ori_rate);
							AddBurn(burn_count,h(load_results.at(1)),stold(load_results.at(2)),stold(load_results.at(3)),stold(load_results.at(4)),burn_ori,burn_ori_rate);
							burn_count++;
						}
						if(line.compare(1, 8,"Addition")==0){
							cout<<"This functionality has not yet been added";
							throw exception();
						}
					}
				}
				file.close();
			}
        }

		void calculate_a_a_dot(vector <long double> *a_local, vector <long double> *a_dot_local,
		// string *body_id_local, vector<long double> *r, vector<long double> *v){
		long double *body_id_local, vector<long double> *r, vector<long double> *v){	
			// #pragma omp parallel for
			for (int i=0; i<body_ids.size(); i++){

				vector <long double> acting_values_local;
				// string acting_id_local;
				long double acting_id_local;
				vector <long double> v_rel_local, r_rel_local;
				long double abs_r_rel_local, r_dot_v_relative_local;
				long double a_comp_local, a_dot_local_comp_1, a_dot_local_comp_2;

				if(body_ids[i] != *body_id_local){
					acting_values_local = body_info[i];
					r_rel_local.clear();
					v_rel_local.clear();
					for( int i=0;i<3;i++){
						r_rel_local.push_back(r->at(i)-acting_values_local.at(1+i));
					}
					for( int i=0;i<3;i++){
						v_rel_local.push_back(v->at(i)-acting_values_local.at(4+i));
					}
					for( int i=0;i<3;i++){
						r_dot_v_relative_local=r_rel_local.at(i)*v_rel_local.at(i);
					}

					abs_r_rel_local = sqrt(pow(r_rel_local.at(0),2)+pow(r_rel_local.at(1),2)+pow(r_rel_local.at(2),2));

					a_comp_local = -G*acting_values_local.at(0)/pow(abs_r_rel_local,3);
					a_dot_local_comp_1 = 3*G*acting_values_local.at(0)*r_dot_v_relative_local/pow(abs_r_rel_local,5);
					a_dot_local_comp_2 = -G*acting_values_local.at(0)/pow(abs_r_rel_local,3);

					for( int i=0;i<3;i++){
						a_local->at(i)=a_local->at(i)+a_comp_local*r_rel_local.at(i);
						a_dot_local->at(i)=a_dot_local->at(i)+a_dot_local_comp_1*r_rel_local.at(i)+a_dot_local_comp_2*v_rel_local.at(i);
					}
				}
			}
		}

		// void CalcA_Burn(vector <long double> *a_local, string *body_id_local){
		void CalcA_Burn(vector <long double> *a_local, long double *body_id_local){
			// omp_set_num_threads(4);
			// #pragma omp parallel for
			// for (pair<string, pair<string, vector<long double> > > burn_itt : burns){
			for(int i=0; i<burn_count; i++){
				// string burn_id_local, burn_body_local;
				long double burn_body_local;
				// pair<string, vector<long double> > burn_values_local;
				pair<long double, vector<long double> > burn_values_local;
				vector <long double> burn_vector_local;
				burn_values_local = burns[i];
				burn_body_local = burn_values_local.first;
				if(burn_body_local == *body_id_local){
					burn_vector_local = burn_values_local.second;
					if(time>burn_vector_local.at(0)&&time<burn_vector_local.at(1)){
						for( int i=0;i<3;i++){
							a_local->at(i) = burn_vector_local.at(2)*burn_vector_local.at(3+i);
						}
					/*
							burn_vector_local_next.clear();
							for( int i=0;i<3;i++){
								burn_vector_local_next.push_back(burn_vector_local.at(i));
							}
							for( int i=0;i<3;i++){
								burn_vector_local_next.push_back(burn_vector_local.at(3+i)+burn_vector_local.at(6+i)*timestep);
							}
							for( int i=0;i<3;i++){
								burn_vector_local_next.push_back(burn_vector_local.at(6+i));
							}
							 * */
					}
				}
			}
		}


		// vector<long double> CalcAAndDot(string *body_id_local, vector<long double> *r, vector<long double> *v){
		vector<long double> CalcAAndDot(long double *body_id_local, vector<long double> *r, vector<long double> *v){
			vector <long double> a_local, a_dot_local;
			vector <long double> output_vec_local;

			for(int i=0;i<3;i++){
				a_local.push_back(0);
				a_dot_local.push_back(0);
			}

			calculate_a_a_dot(&a_local, &a_dot_local, body_id_local, r, v);

			//TODO: Swap round looking at body id first with looking at time
			CalcA_Burn(&a_local, body_id_local);

			output_vec_local.clear();
			for(int i=0;i<3;i++){
				output_vec_local.push_back(a_local.at(i));
			}
			for(int i=0;i<3;i++){
				output_vec_local.push_back(a_dot_local.at(i));
			}
			return output_vec_local;
		}

		void calculate_body(int j){

			// string body_id_local;
			long double body_id_local;
			vector <long double> body_values_local;
			vector <long double> body_final_local;
			vector <long double> r_0_local, v_0_local, accel_output_local;
			vector <long double> a_0_local, a_0_local_dot;
			vector <long double> r_p_local, v_p_local;
			vector <long double> a_p_local, a_p_local_dot;

			body_id_local = body_ids[j];
			body_values_local = body_info[j];

			for( int i=0;i<3;i++){
				r_0_local.push_back(body_values_local.at(1+i));
				v_0_local.push_back(body_values_local.at(4+i));
			}

			accel_output_local = CalcAAndDot(&body_id_local,&r_0_local,&v_0_local);

			for( int i=0;i<3;i++){
				a_0_local.push_back(accel_output_local.at(i));
				a_0_local_dot.push_back(accel_output_local.at(3+i));
			}
			for( int i=0;i<3;i++){
				r_p_local.push_back(r_0_local.at(i)+v_0_local.at(i)*timestep+0.5*a_0_local.at(i)*pow(timestep,2)+(1/6)*a_0_local_dot.at(i)*pow(timestep,3));
				v_p_local.push_back(v_0_local.at(i)+a_0_local.at(i)*timestep+0.5*a_0_local_dot.at(i)*pow(timestep,2));
			}

			for(itts=0;itts<2;itts++){
				accel_output_local = CalcAAndDot(&body_id_local,&r_p_local,&v_p_local);
				a_p_local.clear();a_p_local_dot.clear();
				for( int i=0;i<3;i++){
					a_p_local.push_back(accel_output_local.at(i));
					a_p_local_dot.push_back(accel_output_local.at(3+i));
				}

				for( int i=0;i<3;i++){
					v_p_local.push_back(v_0_local.at(i)+0.5*(a_0_local.at(i)+a_p_local.at(i))*timestep+(1/12)*(a_0_local_dot.at(i)-a_p_local_dot.at(i)*pow(timestep,2)));
					r_p_local.push_back(r_0_local.at(i)+0.5*(v_p_local.at(i)+v_0_local.at(i))*timestep+(1/12)*(a_0_local.at(i)-a_p_local.at(i))*pow(timestep,2));
				}
			}

			body_final_local.push_back(body_values_local.at(0));
			for( int i=0;i<3;i++){
				body_final_local.push_back(r_p_local.at(i));
			}
			for( int i=0;i<3;i++){
				body_final_local.push_back(v_p_local.at(i));
			}
			body_info_next[j]=body_final_local;

		}

		void Step() {
			bodies_next.clear();
			next_timestep = 9999999999;
			
			#pragma omp parallel for
			for(int i=0; i<body_ids.size(); i++){
				calculate_body(i);
			}

			body_info = body_info_next;

			time = time + timestep;
		}

		void Output(){
			if(time_counter==output_rate){
				cout<<"#"+to_string(time)+"\n";
				for (int i=0; i<body_ids.size(); i++){

					// body_id = body_ids[i];
					string body_string_id = body_string_ids.at(i);

					body_values = body_info[i];
					// cout<<body_id;
					cout<<body_string_id;
					for( int i=0;i<body_values.size();i++){
						cout<<","+to_string(body_values.at(i));
					}
					cout<<"\n";
				}
				time_counter=0;
			}
			else{
				time_counter++;
			}
		}
};
