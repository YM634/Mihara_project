//
//  sc_dynamics_main_fixMF.cpp
//  Lib_TDDM_ver2
//
//  Created by Yuta Murakami on 2020/06/28.
//  Copyright © 2020 Yuta Murakami. All rights reserved.
//

#include "sc_dynamics_main_fixMF.hpp"
#include <chrono>

int main(int argc, char *argv[]){

    //..................................................
    //parameters for system
    //..................................................
    double T;
    double Jhop = 1.0;
    double mu;
    double U_el;
    double g_elph;
    double w0_ph;
    int L;
    
    //..................................................
    //parameters for caclucation condition
    //..................................................
    string SC_symmetry; //"S" or "D" is available
    string elph_type; //"S" , "D" or "hop" is available
    //calclulation parameters
    double dt;
    int Nt;

    //..................................................
    // Excitation condition
    //..................................................
    vector<double> Ax_t,Ay_t;
    vector<double> Ex_t,Ey_t;
    //..................................................
    // Physical quantities
    //..................................................
    SC_model_basic model;
    //============================================================================
    //   READ INPUT
    //============================================================================
    
      {
          cout<<"======================================================"<<endl;
          cout<<" Getting parameters"<<endl;
          cout<<"======================================================"<<endl;

          find_param(argv[1],"__T=",T);
          find_param(argv[1],"__mu=",mu);
          find_param(argv[1],"__U_el=",U_el);
          find_param(argv[1],"__g_elph=",g_elph);
          find_param(argv[1],"__w0_ph=",w0_ph);
          
          find_param(argv[1],"__L=",L);
          find_param(argv[1],"__Nt=",Nt);
          find_param(argv[1],"__dt=",dt);
          
          find_param(argv[1],"__SC_symmetry=",SC_symmetry);
          find_param(argv[1],"__elph_type=",elph_type);
          
          //..................................................
          //info of excitations
          //..................................................
          
          find_param_tvector(argv[1],"__Ex_t=",Ex_t,Nt);
          find_param_tvector(argv[1],"__Ax_t=",Ax_t,Nt);
          
          find_param_tvector(argv[1],"__Ey_t=",Ey_t,Nt);
          find_param_tvector(argv[1],"__Ay_t=",Ay_t,Nt);
              
      }
    
    //============================================================================
    //   Initialize quantities
    //============================================================================
    {
         cout<<"======================================================"<<endl;
         cout<<" Initialize DATA and System"<<endl;
         cout<<"======================================================"<<endl;
        model.set_system_params(T, mu, U_el, g_elph, w0_ph, L, SC_symmetry, elph_type);
        model.set_calcluation_params(dt, Nt, Ax_t, Ay_t);
    }
    //============================================================================
    //   Main part
    //============================================================================
    {
        cout<<"======================================================"<<endl;
        cout<<" Main:  "<<endl;
        cout<<"======================================================"<<endl;
        cout<<setprecision(17);
        
        auto eq_start=std::chrono::steady_clock::now();
        model.set_equilibrium();
        auto eq_end=std::chrono::steady_clock::now();
        
        if(model.converged){
            {
                std::ofstream initial("initial_state.tsv"); initial<<std::setprecision(17);
                for(int ik=0;ik<model.Nsys;ik++){
                    auto &r=model.rho_el_k_tstp[ik];
                    initial<<2*real(r(0,1))<<" "<<-2*imag(r(0,1))<<" "<<real(r(0,0)-r(1,1))<<"\n";
                }
            }
            auto evolution_start=std::chrono::steady_clock::now();
            for(int it=0 ; it<Nt ; it++){
                if(it%500==0) cout<<"it= "<<it<<endl;
                model.evaluate_observales_fixMF(it);
                model.time_evoluation_tstp_fixMF(it);
                
            }
            model.evaluate_observales_fixMF(Nt);
            auto evolution_end=std::chrono::steady_clock::now();
            std::ofstream timing("kernel_timing.tsv");timing<<std::setprecision(17);
            timing<<"equilibrium_seconds "<<std::chrono::duration<double>(eq_end-eq_start).count()<<"\n";
            timing<<"evolution_seconds "<<std::chrono::duration<double>(evolution_end-evolution_start).count()<<"\n";
        }
    }
    
    //============================================================================
    //   Print observables
    //============================================================================
    if(model.converged){
    {
        char filename[100];
        
        sprintf(filename, "obs_t.dat");
        ofstream fout_obs(filename);
        fout_obs<<setprecision(17);
        fout_obs<<"#0:time"<<" ";
        fout_obs<<"#1:n_up"<<" ";
        fout_obs<<"#2:n_down"<<" ";
        fout_obs<<"#3:Phi_sc.re"<<" ";
        fout_obs<<"#4:Phi_sc.im"<<" ";
        fout_obs<<"#5:X_ph_q0"<<" ";
        fout_obs<<"#6:Etot"<<" ";
        fout_obs<<"#7:Ekin"<<" ";
        fout_obs<<"#8:Eint"<<" ";
        fout_obs<<"#9:E_nx"<<" ";
        fout_obs<<"#10:E_ph"<<" ";

        fout_obs<<endl;
        
        for(int it=0;it<=Nt;it++){
            
            fout_obs<<(double)it*dt<<" ";
            fout_obs<<model.n_up[it]<<" ";
            fout_obs<<model.n_do[it]<<" ";
            fout_obs<<real(model.Phi_sc[it])<<" ";
            fout_obs<<imag(model.Phi_sc[it])<<" ";
            fout_obs<<model.X_ph_q0[it]<<" ";
            fout_obs<<model.E_tot[it]<<" ";
            fout_obs<<model.E_kin[it]<<" ";
            fout_obs<<model.E_int[it]<<" ";
            fout_obs<<model.E_nx[it]<<" ";
            fout_obs<<model.E_ph[it]<<" ";
            fout_obs<<endl;
        }
        fout_obs.close();
        
    }
    {
        char filename[100];
        
        sprintf(filename, "Jcurr_t.dat");
        ofstream fout_obs(filename);
        fout_obs<<setprecision(17);
        fout_obs<<"#0:time"<<" ";
        fout_obs<<"#1:Jx"<<" ";
        fout_obs<<"#2:Jy"<<" ";
        fout_obs<<endl;
        
        for(int it=0;it<=Nt;it++){
            
            fout_obs<<(double)it*dt<<" ";
            fout_obs<<model.J_curr[it](0)<<" ";
            fout_obs<<model.J_curr[it](1)<<" ";
            fout_obs<<endl;
        }
        fout_obs.close();
        
    }
        {
            char filename[100];
            sprintf(filename, "check_energy.dat");
            ofstream fout_obs(filename);
            fout_obs<<setprecision(17);
            
            double diff_en=0.0;
            
            fout_obs<<0.0<<" ";
            fout_obs<<diff_en<<" ";
            fout_obs<<0.0<<endl;
            
            for(int it=1;it<=Nt;it++){
                Vector2d Et_old,Et_new;
                Et_old(0) = Ex_t[it-1];
                Et_old(1) = Ey_t[it-1];
                Et_new(0) = Ex_t[it];
                Et_new(1) = Ey_t[it];
            
                diff_en += 0.5*(Et_old.dot(model.J_curr[it-1]) + Et_new.dot(model.J_curr[it]))*dt;
                fout_obs<<(double)it*dt<<" ";
                fout_obs<<diff_en<<" ";
                fout_obs<<model.E_tot[it] - model.E_tot[0]<<" ";
                fout_obs<<endl;
            }
            fout_obs.close();
        }
        {
            char filename[100];
            sprintf(filename, "gap_eq.dat");
            ofstream fout_obs(filename);
            fout_obs<<setprecision(17);
            fout_obs<<"# 0:kx 1:ky 2:gap_k 3:gap_k-gap0_k 4:gap0_k"<<endl;
            for(int ikx=0 ; ikx<model.L ; ikx++){
                for(int iky=0 ; iky<model.L ; iky++){
                    int ik = ikx + L*iky;
                    fout_obs<<model.k_vec[ik](0)<<" "<<model.k_vec[ik](1)<<" ";
                    fout_obs<<model.gap_k[ik]<<" "<<model.gap_k[ik]-model.gap0_k[ik]
                    <<" "<<model.gap0_k[ik]<<endl;
                }
                fout_obs<<endl;
            }
            
        }
    }
}
