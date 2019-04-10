
import java.io.*;
import java.util.*;

import org.jlab.groot.math.*;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.base.GStyle;

public class HTCC{
	boolean userTimeBased, write_volatile;
	int runNum;
	boolean[] trigger_bits;
	public float EBeam;
        public int e_part_ind, e_sect, e_track_ind;
        public float RFtime, e_mom, e_theta, e_phi, e_vx, e_vy, e_vz, e_ecal_X, e_ecal_Y, e_ecal_Z, e_ecal_E, e_track_chi2, e_vert_time, e_vert_time_RF, e_Q2, e_xB, e_W;
        public float e_HTCC, e_LTCC, e_pcal_e, e_etot_e, e_TOF_X, e_TOF_Y, e_TOF_Z, e_HTCC_X, e_HTCC_Y, e_HTCC_Z;
        public float e_DCR1_X, e_DCR1_Y, e_DCR1_Z, e_DCR2_X, e_DCR2_Y, e_DCR2_Z, e_DCR3_X, e_DCR3_Y, e_DCR3_Z;
        public float e_DCR1_uX, e_DCR1_uY, e_DCR1_uZ, e_DCR2_uX, e_DCR2_uY, e_DCR2_uZ, e_DCR3_uX, e_DCR3_uY, e_DCR3_uZ;
	public float e_DCR2_the, e_DCR2_phi;

	public H2F[] H_e_theta_mom, H_e_phi_mom, H_e_theta_phi, H_e_vz, H_e_sampl, H_e_vtime, H_e_trk_chi2, H_e_HTCC;
	public H2F[] H_e_Ring_theta, H_e_side_phi;
	public H1F[] H_HTCC_adc, H_HTCC_nphe, H_HTCC2_nphe;
	public HTCC(int reqR, float reqEb, boolean reqTimeBased, boolean reqwrite_volatile){
        	runNum = reqR;userTimeBased=reqTimeBased;
		write_volatile = reqwrite_volatile;
		EBeam = 2.2f;
                if(reqEb>0 && reqEb<4)EBeam=2.22f;
                if(reqEb>4 && reqEb<7.1)EBeam=6.42f;
                if(reqEb>7.1 && reqEb<9)EBeam=7.55f;
                if(reqEb>9)EBeam=10.6f;
		trigger_bits = new boolean[32];
		H_e_theta_mom = new H2F[7];
		H_e_phi_mom = new H2F[7];
		H_e_theta_phi = new H2F[7];
		H_e_vz = new H2F[7];
		H_e_sampl = new H2F[7];
		H_e_vtime = new H2F[7];
		H_e_trk_chi2 = new H2F[7];
		H_e_HTCC = new H2F[7];
		H_e_Ring_theta = new H2F[7];
		H_e_side_phi = new H2F[7];
		for(int s=0;s<7;s++){
			H_e_theta_mom[s] = new H2F(String.format("H_e_S%d_theta_mom",s+1),String.format("e #theta vs p S%d",s+1),100,0,EBeam,100,0,40);
			H_e_theta_mom[s].setTitleX("p (GeV)");
			H_e_theta_mom[s].setTitleY("#theta (^o)");
			H_e_phi_mom[s] = new H2F(String.format("H_e_S%d_phi_mom",s+1),String.format("e #phi vs p S%d",s+1),100,0,EBeam,100,-180,180);
			H_e_phi_mom[s].setTitleX("p (GeV)");
			H_e_phi_mom[s].setTitleY("#phi (^o)");
			H_e_theta_phi[s] = new H2F(String.format("H_e_S%d_theta_phi",s+1),String.format("e #theta vs #phi S%d",s+1),100,-180,180,100,0,40);
			H_e_theta_phi[s].setTitleX("#phi (^o)");
			H_e_theta_phi[s].setTitleY("#theta (^o)");
			H_e_vz[s] = new H2F(String.format("H_e_S%d_vz",s+1),String.format("e vz vs #theta S%d",s+1),100,0,40,100,-25,25);
			H_e_vz[s].setTitleX("#theta (^o)");
			H_e_vz[s].setTitleY("vz (cm)");
			H_e_sampl[s] = new H2F(String.format("H_e_S%d_sampl",s+1),String.format("e sampling vs p S%d",s+1),100,0,EBeam,100,0,0.5);
			H_e_sampl[s].setTitleX("p (GeV)");
			H_e_sampl[s].setTitleY("ECAL/p");
			H_e_vtime[s] = new H2F(String.format("H_e_S%d_vtim",s+1),String.format("e t_v vs p S%d",s+1),100,0,EBeam,120,-1.1,1.1);
			H_e_vtime[s].setTitleX("p (GeV)");
			H_e_vtime[s].setTitleY("t vertex (ns)");
			H_e_trk_chi2[s] = new H2F(String.format("H_e_S%d_trkchi2",s+1),String.format("e #chi^2 vs #theta S%d",s+1),100,0,40,100,0,700);
			H_e_trk_chi2[s].setTitleX("#theta (^o)");
			H_e_trk_chi2[s].setTitleY("DC #chi^2");
			H_e_HTCC[s] = new H2F(String.format("H_e_S%d_HTCC",s+1),String.format("e HTCC vs #theta S%d",s+1),100,0,40,100,0,100);
			H_e_HTCC[s].setTitleX("#theta (^o)");
			H_e_HTCC[s].setTitleY("nphe HTCC");
			H_e_Ring_theta[s] = new H2F(String.format("H_e_Ring_theta_S%d",s+1),String.format("e S%d HTCC ring vs DCR1 #theta",s+1),100,0,40,4,0.5,4.5);
			H_e_Ring_theta[s].setTitleX("#theta (^o)");
			H_e_Ring_theta[s].setTitleY("ring");
			H_e_side_phi[s] = new H2F(String.format("H_e_side_phi_S%d",s+1),String.format("e S%d HTCC side vs DCR1 #phi",s+1),100,-30,30,2,-1,1);
			H_e_side_phi[s].setTitleX("#phi (^o)");
			H_e_side_phi[s].setTitleY("side");
		}
		//others
		H_HTCC_adc = new H1F[48];
		H_HTCC_nphe = new H1F[48];
		H_HTCC2_nphe = new H1F[48];
		for(int r=0;r<4;r++){
			for(int side=0;side<2;side++){
				for(int s=0;s<6;s++){
					int counter = r + 4*( side + 2*s );
					String stringSide = "left";
					if(side==1)stringSide = "right";
					String histitle = String.format("HTCC ADC S%d, Ring %d, %s",s+1,r+1,stringSide);
					H_HTCC_adc[counter] = new H1F(String.format("H_HTCC_adc%d",s+1),histitle,100,0,10000);
					histitle = String.format("HTCC NPHE S%d, Ring %d, %s",s+1,r+1,stringSide);
					H_HTCC_nphe[counter] = new H1F(String.format("H_HTCC_nphe%d",s+1),histitle,100,0,50);
					histitle = String.format("HTCC UNMATCHED NPHE S%d, Ring %d, %s",s+1,r+1,stringSide);
					H_HTCC2_nphe[counter] = new H1F(String.format("H_HTCC2_nphe%d",s+1),histitle,100,0,50);
				}
			}
		}
		//for(int s=0;s<48;s++){
		//	H_HTCC_adc[s] = new H1F(String.format("H_HTCC_adc%d",s+1),String.format("H_HTCC_adc%d",s+1),100,0,10000);
		//	H_HTCC_nphe[s] = new H1F(String.format("H_HTCC_nphe%d",s+1),String.format("H_HTCC_nphe%d",s+1),100,0,50);
		//}
	}
        public int makeElectron(DataBank bank){
                for(int k = 0; k < bank.rows(); k++){
                        int pid = bank.getInt("pid", k);
                        byte q = bank.getByte("charge", k);
                        float px = bank.getFloat("px", k);
                        float py = bank.getFloat("py", k);
                        float pz = bank.getFloat("pz", k);
                        int status = bank.getShort("status", k);
                        boolean inDC = (status>=2000 && status<4000);
                        e_mom = (float)Math.sqrt(px*px+py*py+pz*pz);
			e_theta = (float)Math.toDegrees(Math.acos(pz/e_mom));
			e_vz = bank.getFloat("vz", k);
                        if(inDC && pid == 11 && e_mom>EBeam*0.15 && e_theta>6.5 && Math.abs(e_vz)<20){
                                e_phi = (float)Math.toDegrees(Math.atan2(py,px));
                                e_vx = bank.getFloat("vx", k);
                                e_vy = bank.getFloat("vy", k);
                                return k;
                        }
                }
                return -1;
        }    

        public void getElecEBECal(DataBank bank){
                e_ecal_E=0;
		for(int k = 0; k < bank.rows(); k++){
                        int det = bank.getInt("layer", k);
                        short pind = bank.getShort("pindex",k);
                        if(det==1 && pind==e_part_ind){
                                e_ecal_X = bank.getFloat("x",k);
                                e_ecal_Y = bank.getFloat("y",k);
                                e_ecal_Z = bank.getFloat("z",k);
                                e_ecal_E += bank.getFloat("energy",k);
                                e_pcal_e += bank.getFloat("energy",k);
                                e_sect = bank.getByte("sector",k);
                        }
                        if(det==4 && pind==e_part_ind){
                                e_ecal_E += bank.getFloat("energy",k);
                                e_etot_e += bank.getFloat("energy",k);
                        }
                        if(det==7 && pind==e_part_ind){
                                e_ecal_E += bank.getFloat("energy",k);
                                e_etot_e += bank.getFloat("energy",k);
                        }
                }
        }
        public void getElecEBCC(DataBank bank){
                for(int k = 0; k < bank.rows(); k++){
                        short pind = bank.getShort("pindex",k);
                        if(bank.getByte("detector",k)==15 && pind==e_part_ind){
                                //H_e_HTCC_nphe.fill(bank.getShort("nphe",k));
                                //H_e_HTCC_xy.fill(bank.getFloat("x",k) , bank.getFloat("y",k));
                                //e_HTCC = (float)bank.getShort("nphe",k);
                                e_HTCC = (float)bank.getFloat("nphe",k);
                                e_HTCC_X = bank.getFloat("x",k);
                                e_HTCC_Y = bank.getFloat("y",k);
                                e_HTCC_Z = bank.getFloat("z",k);

                        }    
                        if(bank.getByte("detector",k)==16 && pind==e_part_ind){
                                //H_e_LTCC_nphe.fill(bank.getShort("nphe",k));
                                //H_e_LTCC_xy.fill(bank.getFloat("x",k) , bank.getFloat("y",k));
                                //e_LTCC = (float)bank.getShort("nphe",k);
                                e_LTCC = (float)bank.getFloat("nphe",k);
                        }    
                }    
        }
        public void getElecEBTOF(DataBank bank){
                for(int k = 0; k < bank.rows(); k++){
                        short pind = bank.getShort("pindex",k);
                        if(pind==e_part_ind && bank.getFloat("energy",k)>5){
                                //H_e_TOF_xy.fill(bank.getFloat("x",k) , bank.getFloat("y",k));
                                //H_e_TOF_t_path.fill(bank.getFloat("time",k),bank.getFloat("path",k));
                                e_vert_time = bank.getFloat("time",k) - bank.getFloat("path",k)/29.98f;
                                float time = (e_vert_time-RFtime+1.002f)%2.004f;time -= 1.002f;
                                e_vert_time_RF = time;
                                //H_e_vt1.fill(e_vert_time_RF);
                                //H_e_vt2.fill(time2);
                                e_TOF_X = bank.getFloat("x",k);
                                e_TOF_Y = bank.getFloat("y",k);
                                e_TOF_Z = bank.getFloat("z",k);
                        }    
                }    
        }    
        public void fillEBTrack(DataBank bank){
                e_track_ind=-1;
                for(int k = 0; k < bank.rows(); k++){
                        short pind = bank.getShort("pindex",k);
                        if(pind==e_part_ind){
                                e_track_chi2 =  bank.getFloat("chi2",k);
                                e_track_ind = bank.getShort("index",k);
                        }    
                }    
        }    
        public void getTrigTBTrack(DataBank bank){
                 if(e_track_ind>-1 && e_track_ind<bank.rows() ){
                         e_track_chi2 = bank.getFloat("chi2" , e_track_ind);
                         e_sect = bank.getInt("sector", e_track_ind);
                 }
        }
        public void getTBTrack(DataBank bank){ 
                 if(e_track_ind>-1 && e_track_ind<bank.rows()){
                         e_track_chi2 = bank.getFloat("chi2" , e_track_ind);
                         e_sect = bank.getInt("sector", e_track_ind);
			 e_DCR2_X = bank.getFloat("t1_x" , e_track_ind);
			 e_DCR2_Y = bank.getFloat("t1_y" , e_track_ind);
			 e_DCR2_Z = bank.getFloat("t1_z" , e_track_ind);
			 e_DCR2_uX = bank.getFloat("t1_px" , e_track_ind);
			 e_DCR2_uY = bank.getFloat("t1_py" , e_track_ind);
			 e_DCR2_uZ = bank.getFloat("t1_pz" , e_track_ind);
			 Vector3 vDCR2pos = new Vector3(e_DCR2_X,e_DCR2_Y,e_DCR2_Z);
			 vDCR2pos.rotateZ( -3.141597f*(e_sect-1)/3f );
			 e_DCR2_the = (float)Math.toDegrees(vDCR2pos.theta());
			 e_DCR2_phi = (float)Math.toDegrees(vDCR2pos.phi());
                 }    
        }
	public void fill_eHists(){
		H_e_theta_mom[e_sect-1].fill(e_mom,e_theta);
		H_e_phi_mom[e_sect-1].fill(e_mom,e_phi);
		H_e_theta_phi[e_sect-1].fill(e_phi,e_theta);
		H_e_vz[e_sect-1].fill(e_theta,e_vz);
		H_e_sampl[e_sect-1].fill(e_mom,e_ecal_E/e_mom);
		H_e_vtime[e_sect-1].fill(e_mom,e_vert_time_RF);
		H_e_trk_chi2[e_sect-1].fill(e_theta,e_track_chi2);
		H_e_HTCC[e_sect-1].fill(e_theta,e_HTCC);
		H_e_theta_mom[6].fill(e_mom,e_theta);
		H_e_phi_mom[6].fill(e_mom,e_phi);
		H_e_theta_phi[6].fill(e_phi,e_theta);
		H_e_vz[6].fill(e_theta,e_vz);
		H_e_sampl[6].fill(e_mom,e_ecal_E/e_mom);
		H_e_vtime[6].fill(e_mom,e_vert_time_RF);
		H_e_trk_chi2[6].fill(e_theta,e_track_chi2);
		H_e_HTCC[6].fill(e_theta,e_HTCC);
	}
	public void HTCCadc(DataBank bank){
		int bestmatch = -1;
		float matchnadc = 0;
		for(int r=0;r<bank.rows();r++){
			int sect = bank.getInt("sector",r);
			int side = bank.getInt("layer",r);
			int ring = bank.getInt("component",r);
			int adcv = bank.getInt("ADC",r);
			int peds = bank.getInt("ped",r);
			int counter = (ring-1) + 4*( (side-1) + 2*(sect-1) );
			if(counter>47)System.out.println("ERROR COUNTING PMTS "+counter);
			//else H_HTCC_adc[counter].fill(adcv-peds+100);
			else if(e_sect==sect){
				boolean geomMatch = false;
				geomMatch = ( Math.abs( 1f*ring - (1f + 3f*(e_DCR2_the-10f)/22.5f) ) < 0.75f) && (e_DCR2_phi*(side-1.5f)>0);
				if(geomMatch && adcv>matchnadc){
					//H_HTCC_adc[counter].fill(adcv+100);
					matchnadc = adcv;
					bestmatch = r;
				}
			}
		}
		if(bestmatch>-1){
			int r = bestmatch;
			int counter = (bank.getInt("component",r)-1) + 4*( (bank.getInt("layer",r)-1) + 2*(bank.getInt("sector",r)-1) );
			H_HTCC_adc[counter].fill(bank.getInt("ADC",r)+100);
		}
	}
	public void HTCCrec(DataBank bank, DataBank hbank){
		int bestmatch = -1;
		float matchnphe = 0;
		for(int r=0;r<bank.rows();r++){
			int IDhit = bank.getInt("id",r);
			int nhits = bank.getInt("nhits",r);
			float nphe = bank.getFloat("nphe",r);
			if(nhits==1){
				int sect = hbank.getInt("sector",IDhit);
				int side = hbank.getInt("layer",IDhit);
				int ring = hbank.getInt("component",IDhit);
				int counter = (ring-1) + 4*( (side-1) + 2*(sect-1) );
				if(e_sect==sect){
					boolean geomMatch = false;
					geomMatch = ( Math.abs( 1f*ring - (1f + 3f*(e_DCR2_the-10f)/22.5f) ) < 0.33f) && (e_DCR2_phi*(side-1.5f)>0);
					if(geomMatch && nphe>matchnphe){
						matchnphe=nphe;bestmatch=r;
					}
				}
			}
		}
		if(bestmatch>-1){
			int r = bestmatch;
			int IDhit = bank.getInt("id",r);
			int nhits = bank.getInt("nhits",r);
			float nphe = bank.getFloat("nphe",r);
			int sect = hbank.getInt("sector",IDhit);
			int side = hbank.getInt("layer",IDhit);
			int ring = hbank.getInt("component",IDhit);
			int counter = (ring-1) + 4*( (side-1) + 2*(sect-1) );
			H_HTCC_nphe[counter].fill(nphe);
			H_e_Ring_theta[sect-1].fill(e_DCR2_the,ring);
			H_e_side_phi[sect-1].fill(e_DCR2_phi,side-1.5f);
			H_e_Ring_theta[6].fill(e_DCR2_the,ring);
			H_e_side_phi[6].fill(e_DCR2_phi,side-1.5f);
		
		}
		for(int r=0;r<bank.rows();r++){
			if(r!=bestmatch){
				int IDhit = bank.getInt("id",r);
				float nphe = bank.getFloat("nphe",r);
				int sect = hbank.getInt("sector",IDhit);
				int side = hbank.getInt("layer",IDhit);
				int ring = hbank.getInt("component",IDhit);
				int counter = (ring-1) + 4*( (side-1) + 2*(sect-1) );
				H_HTCC2_nphe[counter].fill(nphe);
			}
		}
	}
	public void processEvent(DataEvent event){
		e_part_ind = -1;
		RFtime=0;
		if(event.hasBank("RUN::config")){
			DataBank confbank = event.getBank("RUN::config");
			long TriggerWord = confbank.getLong("trigger",0);
			for (int i = 31; i >= 0; i--) {trigger_bits[i] = (TriggerWord & (1 << i)) != 0;} 
			if(event.hasBank("RUN::rf")){
				for(int r=0;r<event.getBank("RUN::rf").rows();r++){
					if(event.getBank("RUN::rf").getInt("id",r)==1)RFtime=event.getBank("RUN::rf").getFloat("time",r);
				}    
			}
                DataBank partBank = null, trackBank = null, trackDetBank = null, ecalBank = null, cherenkovBank = null, scintillBank = null;
                if(userTimeBased){
                        if(event.hasBank("REC::Particle"))partBank = event.getBank("REC::Particle");
                        if(event.hasBank("REC::Track"))trackBank = event.getBank("REC::Track");
                        if(event.hasBank("TimeBasedTrkg::TBTracks"))trackDetBank = event.getBank("TimeBasedTrkg::TBTracks");
                        if(event.hasBank("REC::Calorimeter")) ecalBank = event.getBank("REC::Calorimeter");
                        if(event.hasBank("REC::Cherenkov"))cherenkovBank = event.getBank("REC::Cherenkov");
                        if(event.hasBank("REC::Scintillator"))scintillBank = event.getBank("REC::Scintillator");
                }
                if(!userTimeBased){
                        if(event.hasBank("RECHB::Particle"))partBank = event.getBank("RECHB::Particle");
                        if(event.hasBank("RECHB::Track"))trackBank = event.getBank("RECHB::Track");
                        if(event.hasBank("HitBasedTrkg::HBTracks"))trackDetBank = event.getBank("HitBasedTrkg::HBTracks");
                        if(event.hasBank("RECHB::Calorimeter")) ecalBank = event.getBank("RECHB::Calorimeter");
                        if(event.hasBank("RECHB::Cherenkov"))cherenkovBank = event.getBank("RECHB::Cherenkov");
                        if(event.hasBank("RECHB::Scintillator"))scintillBank = event.getBank("RECHB::Scintillator");
                }

			if( (trigger_bits[1] || trigger_bits[2] || trigger_bits[3] || trigger_bits[4] || trigger_bits[5] || trigger_bits[6]) && partBank!=null)e_part_ind = makeElectron(partBank);
			if(e_part_ind==-1)return;
			if(trackBank!=null)fillEBTrack(trackBank);
			if(ecalBank!=null)getElecEBECal(ecalBank);
			if(cherenkovBank!=null)getElecEBCC(cherenkovBank);
			if(scintillBank!=null)getElecEBTOF(scintillBank);
			if(trackDetBank!=null)getTBTrack(trackDetBank);
			float Q2 = (float) (4 * EBeam * e_mom * Math.pow(Math.sin(Math.toRadians(e_theta) / 2),2));
			if(e_mom>EBeam*0.15 && e_theta>6.5 && Q2>0.2 *EBeam/11 && Math.abs(e_vz)<20 && e_track_chi2<500 && e_ecal_E/e_mom>0.15){
				fill_eHists();
				if(event.hasBank("HTCC::adc"))HTCCadc(event.getBank("HTCC::adc"));
				if(event.hasBank("HTCC::rec")&&event.hasBank("HTCC::adc"))HTCCrec(event.getBank("HTCC::rec"),event.getBank("HTCC::adc"));
			}
		}
	}
        public void plot() {
		EmbeddedCanvas can_e_HTCC  = new EmbeddedCanvas();
		can_e_HTCC.setSize(3500,4000);
		can_e_HTCC.divide(7,8);
		can_e_HTCC.setAxisTitleSize(18);
		can_e_HTCC.setAxisFontSize(18);
		can_e_HTCC.setTitleSize(18);
		for(int s=0;s<7;s++){
			can_e_HTCC.cd(s);can_e_HTCC.draw(H_e_theta_mom[s]);
			can_e_HTCC.cd(7+s);can_e_HTCC.draw(H_e_phi_mom[s]);
			can_e_HTCC.cd(14+s);can_e_HTCC.draw(H_e_theta_phi[s]);
			can_e_HTCC.cd(21+s);can_e_HTCC.draw(H_e_vz[s]);
			can_e_HTCC.cd(28+s);can_e_HTCC.draw(H_e_sampl[s]);
			can_e_HTCC.cd(35+s);can_e_HTCC.draw(H_e_vtime[s]);
			can_e_HTCC.cd(42+s);can_e_HTCC.draw(H_e_trk_chi2[s]);
			can_e_HTCC.cd(49+s);can_e_HTCC.draw(H_e_HTCC[s]);
		}
		if(runNum>0){
			if(!write_volatile)can_e_HTCC.save(String.format("plots"+runNum+"/HTCC_e.png"));
			if(write_volatile)can_e_HTCC.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/HTCC_e.png"));
			System.out.println(String.format("saved plots"+runNum+"/HTCC_e.png"));
		}
		else{
			can_e_HTCC.save(String.format("plots/HTCC_e.png"));
			System.out.println(String.format("saved plots/HTCC_e.png"));
		}

                EmbeddedCanvas can_HTCC_adc  = new EmbeddedCanvas();
                can_HTCC_adc.setSize(4000,3000);
                can_HTCC_adc.divide(8,6);
                can_HTCC_adc.setAxisTitleSize(18);
                can_HTCC_adc.setAxisFontSize(18);
                can_HTCC_adc.setTitleSize(18);
		for(int s=0;s<48;s++){
			can_HTCC_adc.cd(s);can_HTCC_adc.draw(H_HTCC_adc[s]);
		}
                if(runNum>0){
			if(!write_volatile)can_HTCC_adc.save(String.format("plots"+runNum+"/HTCC_adc.png"));
			if(write_volatile)can_HTCC_adc.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/HTCC_adc.png"));
			System.out.println(String.format("saved plots"+runNum+"/HTCC_adc.png"));
		}
		else{
			can_HTCC_adc.save(String.format("plots/HTCC_adc.png"));
			System.out.println(String.format("saved plots/HTCC_adc.png"));
		}
                EmbeddedCanvas can_HTCC_nphe  = new EmbeddedCanvas();
                can_HTCC_nphe.setSize(4000,3000);
                can_HTCC_nphe.divide(8,6);
                can_HTCC_nphe.setAxisTitleSize(18);
                can_HTCC_nphe.setAxisFontSize(18);
                can_HTCC_nphe.setTitleSize(18);
		for(int s=0;s<48;s++){
			can_HTCC_nphe.cd(s);can_HTCC_nphe.draw(H_HTCC_nphe[s]);
		}
                if(runNum>0){
			if(!write_volatile)can_HTCC_nphe.save(String.format("plots"+runNum+"/HTCC_nphe.png"));
			if(write_volatile)can_HTCC_nphe.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/HTCC_nphe.png"));
			System.out.println(String.format("saved plots"+runNum+"/HTCC_nphe.png"));
		}
		else{
			can_HTCC_nphe.save(String.format("plots/HTCC_nphe.png"));
			System.out.println(String.format("saved plots/HTCC_nphe.png"));
		}
                EmbeddedCanvas can_HTCC2_nphe  = new EmbeddedCanvas();
                can_HTCC2_nphe.setSize(4000,3000);
                can_HTCC2_nphe.divide(8,6);
                can_HTCC2_nphe.setAxisTitleSize(18);
                can_HTCC2_nphe.setAxisFontSize(18);
                can_HTCC2_nphe.setTitleSize(18);
		for(int s=0;s<48;s++){
			can_HTCC2_nphe.cd(s);can_HTCC2_nphe.draw(H_HTCC2_nphe[s]);
		}
                if(runNum>0){
			if(!write_volatile)can_HTCC2_nphe.save(String.format("plots"+runNum+"/HTCC2_nphe.png"));
			if(write_volatile)can_HTCC2_nphe.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/HTCC2_nphe.png"));
			System.out.println(String.format("saved plots"+runNum+"/HTCC2_nphe.png"));
		}
		else{
			can_HTCC2_nphe.save(String.format("plots/HTCC2_nphe.png"));
			System.out.println(String.format("saved plots/HTCC2_nphe.png"));
		}
	}
	public static void main(String[] args) {
                System.setProperty("java.awt.headless", "true");
                GStyle.setPalette("kRainBow");
                int count = 0;
                int runNum = 0;
		boolean useTB = true;
		boolean useVolatile = false;
                String filelist = "list_of_files.txt";
                if(args.length>0)runNum=Integer.parseInt(args[0]);
                if(args.length>1)filelist = args[1];
                int maxevents = 20000000;
                if(args.length>2)maxevents=Integer.parseInt(args[2]);
		float Eb =6.42f;//10.6f;
                if(args.length>3)Eb=Float.parseFloat(args[3]);
		if(args.length>4)if(Integer.parseInt(args[4])==0)useTB=false;
                HTCC ana = new HTCC(runNum,Eb,useTB,useVolatile);
                List<String> toProcessFileNames = new ArrayList<String>();
                File file = new File(filelist);
                Scanner read;
                try {
                        read = new Scanner(file);
                        do { 
                                String filename = read.next();
                                toProcessFileNames.add(filename);

                        }while (read.hasNext());
                        read.close();
                }catch(IOException e){ 
                        e.printStackTrace();
                }   
                int progresscount=0;int filetot = toProcessFileNames.size();
                for (String runstrg : toProcessFileNames) if( count<maxevents ){
                        progresscount++;
                        System.out.println(String.format(">>>>>>>>>>>>>>>> %s",runstrg));
                        File varTmpDir = new File(runstrg);
                        if(!varTmpDir.exists()){System.out.println("FILE DOES NOT EXIST");continue;}
                        System.out.println("READING NOW "+runstrg);
                        HipoDataSource reader = new HipoDataSource();
                        reader.open(runstrg);
                        int filecount = 0;
                        while(reader.hasEvent()&& count<maxevents ) { 
                                DataEvent event = reader.getNextEvent();
                                ana.processEvent(event);
                                filecount++;count++;
                                if(count%10000 == 0) System.out.println(count/1000 + "k events (this is HTCC analysis on "+runstrg+") ; progress : "+progresscount+"/"+filetot);
                        }   
                        reader.close();
                }   
                System.out.println("Total : " + count + " events");
                ana.plot();
        }   
}
