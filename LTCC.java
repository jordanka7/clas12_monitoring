
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

public class LTCC{
	boolean userTimeBased, write_volatile;
	int runNum;
	boolean[] trigger_bits;
	public float EBeam;
        public int e_part_ind, e_sect, e_track_ind, pip_part_ind, pipm_part_ind, pip_sect, pim_sect;
        public float RFtime, e_mom, e_theta, e_phi, e_vx, e_vy, e_vz, e_ecal_X, e_ecal_Y, e_ecal_Z, e_ecal_E, e_track_chi2, e_vert_time, e_vert_time_RF, e_Q2, e_xB, e_W;
        public float e_HTCC, e_LTCC, e_pcal_e, e_etot_e, e_TOF_X, e_TOF_Y, e_TOF_Z, e_HTCC_X, e_HTCC_Y, e_HTCC_Z;
	public float pip_mom, pip_theta, pip_phi, pip_vx, pip_vy, pip_vz, pip_LTCC_X, pip_LTCC_Y, pip_LTCC_Z;
        public float pim_mom, pim_theta, pim_phi, pim_vx, pim_vy, pim_vz, pim_LTCC_X, pim_LTCC_Y, pim_LTCC_Z;
        public float e_DCR1_X, e_DCR1_Y, e_DCR1_Z, e_DCR2_X, e_DCR2_Y, e_DCR2_Z, e_DCR3_X, e_DCR3_Y, e_DCR3_Z;
        public float e_DCR1_uX, e_DCR1_uY, e_DCR1_uZ, e_DCR2_uX, e_DCR2_uY, e_DCR2_uZ, e_DCR3_uX, e_DCR3_uY, e_DCR3_uZ;
	public float e_DCR2_the, e_DCR2_phi;

	public H2F[] H_e_theta_mom, H_e_phi_mom, H_e_theta_phi, H_e_vz, H_e_sampl, H_e_vtime, H_e_trk_chi2, H_e_LTCC;
	public H2F[] H_e_Ring_theta, H_e_side_phi;
	public H1F[] H_LTCC_adc, H_LTCC_nphe, H_LTCC2_nphe;
	public H2F[] H_pion_nphePMT, H_e_nphePMT; 
	public H2F H_e_nphe0_LTCC_XY, H_e_nphe2_LTCC_XY, H_pi_nphe0_LTCC_XY, H_pi_nphe2_LTCC_XY;
	public H1F[] H_LTCC_PMTocc;

	public LTCC(int reqR, float reqEb, boolean reqTimeBased, boolean reqwrite_volatile){
        	runNum = reqR;userTimeBased=reqTimeBased;
		write_volatile = reqwrite_volatile;
		EBeam = 2.2f;
                if(reqEb>0 && reqEb<4)EBeam=2.22f;
                if(reqEb>4 && reqEb<7.1)EBeam=6.42f;
                if(reqEb>7.1 && reqEb<9)EBeam=7.55f;
                if(reqEb>9)EBeam=10.6f;
		trigger_bits = new boolean[32];
		H_pion_nphePMT = new H2F[6];
		H_e_nphePMT = new H2F[6];
		H_LTCC_PMTocc = new H1F[6];

		for(int s=0;s<6;s++){
			H_pion_nphePMT[s] = new H2F(String.format("H_pion_S%d_nphe_vs_PMT",s+1),String.format("Nphe vs PMT (pions) S%d",s+1),36,0.5,36.5,21,-0.5,20.5);
                        H_pion_nphePMT[s].setTitleX("PMT No");
                        H_pion_nphePMT[s].setTitleY("Nphe");
			H_e_nphePMT[s] = new H2F(String.format("H_electron_S%d_nphe_vs_PMT",s+1),String.format("Nphe vs PMT (electrons) S%d",s+1),36,0.5,36.5,21,-0.5,20.5);
                        H_e_nphePMT[s].setTitleX("PMT No");
                        H_e_nphePMT[s].setTitleY("Nphe");
			H_LTCC_PMTocc[s] = new H1F(String.format("H_LTCC_PMT_occ_S%d",s+1),String.format("LTCC PMT Raw Number of Hits (S%d)",s+1),36,0.5,36.5);
                	H_LTCC_PMTocc[s].setTitleX("PMT");
                	H_LTCC_PMTocc[s].setTitleY("events");
		}
		H_e_nphe0_LTCC_XY = new H2F("H_e_nphe0_LTCC_XY","LTCC XY Electrons (nphe>0)",100,-400.,400.,100,-400.,400.);
		H_e_nphe2_LTCC_XY = new H2F("H_e_nphe2_LTCC_XY","LTCC XY Electrons (nphe>2)",100,-400.,400.,100,-400.,400.);
		H_pi_nphe0_LTCC_XY = new H2F("H_e_nphe0_LTCC_XY","LTCC XY Pions (nphe>0)",100,-400.,400.,100,-400.,400.);
                H_pi_nphe2_LTCC_XY = new H2F("H_e_nphe2_LTCC_XY","LTCC XY Pions (nphe>2)",100,-400.,400.,100,-400.,400.);
	}

	public void fill_PMT_Hits(DataBank bankadc) {
		for(int k = 0; k < bankadc.rows(); k++){
			int pmt=0;
                        int sect = bankadc.getInt("sector",k)-1;
			if (bankadc.getByte("layer",k)!=1)System.out.println("Layer = "+bankadc.getByte("layer",k)+" Order = "+bankadc.getByte("order",k)+" Sector = "+sect);
			if (bankadc.getByte("order",k) == 0) pmt = bankadc.getShort("component",k);
			if (bankadc.getByte("order",k) == 1) pmt = bankadc.getShort("component",k)+18;
			H_LTCC_PMTocc[sect].fill(pmt);
		}
	}


        public int isLTCCmatch(DataBank LTCCbank, int index){
		int indexltcc=-1;
		if(userTimeBased){
			//System.out.println("REC Cherenkov N of rows = "+LTCCbank.rows());
			for(int l = 0; l < LTCCbank.rows(); l++) {
				//System.out.println("LTCC row = "+l+" LTCC pindex = "+LTCCbank.getShort("pindex",l)+"Input index from particle bank = "+index+" LTCC Detector = "+LTCCbank.getByte("detector",l));
                                if(LTCCbank.getShort("pindex",l)==index && LTCCbank.getByte("detector",l)==16){
					//System.out.println("LTCC row = "+l+" LTCC pindex = "+LTCCbank.getShort("pindex",l)+" LTCC Detector = "+LTCCbank.getByte("detector",l));
                                        indexltcc=l;
					//System.out.println("Returned value of LTCC index = "+indexltcc);
                                }
                         }
                }
                return indexltcc;
        }

	public int isDCmatch(DataEvent event, int index){
        	int sectordc=-1;
                if(userTimeBased && event.hasBank("REC::Track")){
                        DataBank DCbank = event.getBank("REC::Track");
                                for(int l = 0; l < DCbank.rows(); l++) {
                                        if(DCbank.getShort("pindex",l)==index && DCbank.getInt("detector",l)==6){
                                                        sectordc=DCbank.getByte("sector",l);
                                        }
                                }
                }
                return sectordc;
        }

	public void fill_pion_PMT_Histos(DataEvent event, DataBank bank, DataBank bankrecdet, DataBank bankadc, int id) {
		int det =0;
		int sector = 0;
		int [] sect;
		float nphe = 0;
		int phi = -100;
		int theta = -100;
		int pid = -100;
		int [] countpion;
		int [] counttracks;
		countpion = new int[6];
		counttracks = new int[6];
		int charge = -100;
		int [] pion_index, index;
		pion_index = new int[6];
		index = new int[6];
		sect = new int[6];
		for (int j=0;j<bank.rows();j++) {
			int status = bank.getShort("status", j);
			pid = bank.getInt("pid", j);
			charge = bank.getByte("charge", j);;
			//System.out.println(j+" Particle status = "+status+" pid = "+pid+" charge = "+charge+ "LTCC index = "+isLTCCmatch(bankrecdet,j));
			if ((status>=2000 && status<4000) && isLTCCmatch(bankrecdet,j) != -1 && (pid == id)) {
				//System.out.println(j+" Particle status = "+status+" pid = "+pid+" charge = "+charge+ "LTCC index = "+isLTCCmatch(bankrecdet,j));
				//sector = bankrecdet.getByte("sector",isLTCCmatch(bankrecdet,j));
				sector = isDCmatch(event,j);
				//System.out.println(j+" Particle status = "+status+" pid = "+pid+" charge = "+charge+ "LTCC index = "+isLTCCmatch(bankrecdet,j)+" Sector = "+isDCmatch(event,j)+" "+sector);
				//System.out.println("Sector = "+isDCmatch(event,j)+" "+sector);
				if (sector!=0) {
					countpion[sector-1]++;			
					//System.out.println(sector+" "+countpion[sector-1]);
					pion_index[sector-1] = j; 
					index[sector-1] = isLTCCmatch(bankrecdet,j);
					sect[sector-1] = sector;
				}
			} 
			if ((status>=2000 && status<4000) && isLTCCmatch(bankrecdet, j) != -1 && charge !=0) {
				sector = isDCmatch(event,j);
				if (sector!=0) {
					//System.out.println(+j+" "+sector);
					counttracks[sector-1]++;
				}
			}
			
		}	
		for (int i=0;i<6;i++){
			int pmt = 0;
			if (counttracks[i] == 1 && countpion[i] == 1) {
			//System.out.println("Sector = "+(i+1)+" N of pions = "+countpion[i]+" N of charged tracks = "+counttracks[i]);
				for (int j=0;j<bankadc.rows();j++) {
					if (bankadc.getInt("sector",j) == sect[i]) {
						if (bankadc.getByte("order",j) == 0) pmt = bankadc.getShort("component",j);
						if (bankadc.getByte("order",j) == 1) pmt = bankadc.getShort("component",j)+18;
						nphe = bankrecdet.getFloat("nphe",index[i]);
						H_pion_nphePMT[i].fill(pmt,nphe);
						//System.out.println("sector = "+(i+1)+" adc bank row = "+j+" pmt = "+pmt+" Nphe = "+nphe);
					}	
				}
			}

		}
	}

        public void fill_electron_PMT_Histos(DataEvent event, DataBank bank, DataBank bankrecdet, DataBank bankadc, int id) {
                int det =0;
                int sector = 0;
                int [] sect;
                float nphe = 0;
                int phi = -100;
                int theta = -100;
                int pid = -100;
                int [] countelectron;
                int [] counttracks;
                countelectron = new int[6];
                counttracks = new int[6];
                int charge = -100;
                int [] electron_index, index;
                electron_index = new int[6];
                index = new int[6];
                sect = new int[6];
                for (int j=0;j<bank.rows();j++) {
                        int status = bank.getShort("status", j);
                        pid = bank.getInt("pid", j);
                        charge = bank.getByte("charge", j);
                        //System.out.println(j+" Particle status = "+status+" pid = "+pid+" charge = "+charge+ "LTCC index = "+isLTCCmatch(bankrecdet,j));
                        if ((status>=2000 && status<4000) && isLTCCmatch(bankrecdet,j) != -1 && (pid == id)) {
                                //System.out.println(j+" Particle status = "+status+" pid = "+pid+" charge = "+charge+ "LTCC index = "+isLTCCmatch(bankrecdet,j));
                                //sector = bankrecdet.getByte("sector",isLTCCmatch(bankrecdet,j));
                                sector = isDCmatch(event,j);
                                //System.out.println(j+" Particle status = "+status+" pid = "+pid+" charge = "+charge+ "LTCC index = "+isLTCCmatch(bankrecdet,j)+" Sector = "+isDCmatch(event,j)+" "+sector);
                                //System.out.println("Sector = "+isDCmatch(event,j)+" "+sector);
                                if (sector!=0) {
                                        countelectron[sector-1]++;
                                        //System.out.println(sector+" "+countpion[sector-1]);
                                        electron_index[sector-1] = j;
                                        index[sector-1] = isLTCCmatch(bankrecdet,j);
                                        sect[sector-1] = sector;
                                }
                        }
                        if ((status>=2000 && status<4000) && isLTCCmatch(bankrecdet, j) != -1 && charge !=0) {
                                sector = isDCmatch(event,j);
                                if (sector!=0) {
                                        //System.out.println(+j+" "+sector);
                                        counttracks[sector-1]++;
                                }
                        }

                }
                for (int i=0;i<6;i++){
                        int pmt = 0;
                        if (counttracks[i] == 1 && countelectron[i] == 1) {
                        //System.out.println("Sector = "+(i+1)+" N of pions = "+countpion[i]+" N of charged tracks = "+counttracks[i]);
				for (int j=0;j<bankadc.rows();j++) {
					if (bankadc.getInt("sector",j) == sect[i]) {
						if (bankadc.getByte("order",j) == 0) pmt = bankadc.getShort("component",j);
						if (bankadc.getByte("order",j) == 1) pmt = bankadc.getShort("component",j)+18;
						nphe = bankrecdet.getFloat("nphe",index[i]);
						H_e_nphePMT[i].fill(pmt,nphe);
						//System.out.println("sector = "+(i+1)+" adc bank row = "+j+" pmt = "+pmt+" Nphe = "+nphe);
                                        }
                                }
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
		DataBank trajBank = null, ltccadcBank = null, ltccClusters = null;
                if(userTimeBased){
                        if(event.hasBank("REC::Particle"))partBank = event.getBank("REC::Particle");
                        if(event.hasBank("REC::Track"))trackBank = event.getBank("REC::Track");
                        if(event.hasBank("TimeBasedTrkg::TBTracks"))trackDetBank = event.getBank("TimeBasedTrkg::TBTracks");
                        if(event.hasBank("REC::Calorimeter")) ecalBank = event.getBank("REC::Calorimeter");
                        if(event.hasBank("REC::Cherenkov"))cherenkovBank = event.getBank("REC::Cherenkov");
                        if(event.hasBank("REC::Scintillator"))scintillBank = event.getBank("REC::Scintillator");
			if(event.hasBank("REC::Traj"))trajBank = event.getBank("REC::Traj");
			if(event.hasBank("LTCC::adc"))ltccadcBank = event.getBank("LTCC::adc");
			if(event.hasBank("LTCC::clusters"))ltccClusters = event.getBank("LTCC::clusters");
                }
                if(!userTimeBased){
                        if(event.hasBank("RECHB::Particle"))partBank = event.getBank("RECHB::Particle");
                        if(event.hasBank("RECHB::Track"))trackBank = event.getBank("RECHB::Track");
                        if(event.hasBank("HitBasedTrkg::HBTracks"))trackDetBank = event.getBank("HitBasedTrkg::HBTracks");
                        if(event.hasBank("RECHB::Calorimeter")) ecalBank = event.getBank("RECHB::Calorimeter");
                        if(event.hasBank("RECHB::Cherenkov"))cherenkovBank = event.getBank("RECHB::Cherenkov");
                        if(event.hasBank("RECHB::Scintillator"))scintillBank = event.getBank("RECHB::Scintillator");
			if(event.hasBank("RECHB::Traj"))trajBank = event.getBank("RECHB::Traj");
			if(event.hasBank("LTCC::adc"))ltccadcBank = event.getBank("LTCC::adc");
			if(event.hasBank("LTCC::clusters"))ltccClusters = event.getBank("LTCC::clusters");
                }

			//if (ltccClusters!=null)System.out.println("LTCC Clusters Bank exists");
			if(partBank!=null && cherenkovBank!=null && ltccadcBank!=null){
				fill_pion_PMT_Histos(event,partBank,cherenkovBank,ltccadcBank,211);
				fill_pion_PMT_Histos(event,partBank,cherenkovBank,ltccadcBank,-211);
			}
                        if(partBank!=null && cherenkovBank!=null && ltccadcBank!=null){
                                fill_electron_PMT_Histos(event,partBank,cherenkovBank,ltccadcBank,11);
                                fill_electron_PMT_Histos(event,partBank,cherenkovBank,ltccadcBank,-11);
                        }
			if(ltccadcBank!=null)fill_PMT_Hits(ltccadcBank);
			if(trajBank!=null && partBank!=null && cherenkovBank!=null) fillTraj_LTCC(trajBank,partBank,cherenkovBank);
		}
	}

        public void fillTraj_LTCC(DataBank trajBank, DataBank part, DataBank ltcc){
		for(int i=0;i<part.rows();i++) {
                	float nphe = 0;
                	int pid = -100;
                	int charge = -100;

                        int status = part.getShort("status", i);
                        pid = part.getInt("pid", i);
                        charge = part.getByte("charge", i);;
                        if ((status>=2000 && status<4000) && isLTCCmatch(ltcc,i) != -1 && (pid == 11 || pid == -11 || pid == 211 || pid == -211)) {
				nphe = ltcc.getFloat("nphe",isLTCCmatch(ltcc,i));
				//System.out.println("Nphe = "+nphe+" Traj nrows = "+trajBank.rows());
				for(int r=0;r<trajBank.rows();r++){
                		        if(trajBank.getShort("pindex",r)==i){
                                		if(trajBank.getShort("detId",r) == 43) {
                                        		float e_LTCC_tX = trajBank.getFloat("x",r);
                                        		float e_LTCC_tY = trajBank.getFloat("y",r);
                                        		float e_LTCC_tZ = trajBank.getFloat("z",r);
							if ((pid == 211 || pid == -211) && nphe > 0) H_pi_nphe0_LTCC_XY.fill(e_LTCC_tX,e_LTCC_tY);
							if ((pid == 211 || pid == -211) && nphe > 2) H_pi_nphe2_LTCC_XY.fill(e_LTCC_tX,e_LTCC_tY);
							if ((pid == 11 || pid == -11) && nphe > 0) H_e_nphe0_LTCC_XY.fill(e_LTCC_tX,e_LTCC_tY);
							if ((pid == 11 || pid == -11) && nphe > 2) H_e_nphe2_LTCC_XY.fill(e_LTCC_tX,e_LTCC_tY);
                                		}
                        		}
                		}
			}
			
		} 
        }


        public void plot() {
		EmbeddedCanvas can_e_LTCC  = new EmbeddedCanvas();
		can_e_LTCC.setSize(3500,4000);
		can_e_LTCC.divide(6,10);
		can_e_LTCC.setAxisTitleSize(18);
		can_e_LTCC.setAxisFontSize(18);
		can_e_LTCC.setTitleSize(18);
		for(int s=0;s<6;s++){
			can_e_LTCC.cd(s);can_e_LTCC.draw(H_pion_nphePMT[s]);
			can_e_LTCC.cd(s+6);can_e_LTCC.draw(H_e_nphePMT[s]);
			can_e_LTCC.cd(s+12);can_e_LTCC.draw(H_LTCC_PMTocc[s]);
		}
		can_e_LTCC.cd(18);can_e_LTCC.draw(H_pi_nphe0_LTCC_XY);
		can_e_LTCC.cd(19);can_e_LTCC.draw(H_pi_nphe2_LTCC_XY);
		can_e_LTCC.cd(20);can_e_LTCC.draw(H_e_nphe0_LTCC_XY);
		can_e_LTCC.cd(21);can_e_LTCC.draw(H_e_nphe2_LTCC_XY);
		if(runNum>0){
			if(!write_volatile)can_e_LTCC.save(String.format("plots"+runNum+"/LTCC.png"));
			if(write_volatile)can_e_LTCC.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/LTCC.png"));
			System.out.println(String.format("saved plots"+runNum+"/LTCC.png"));
		}
		else{
			can_e_LTCC.save(String.format("plots/LTCC.png"));
			System.out.println(String.format("saved plots/LTCC.png"));
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
                LTCC ana = new LTCC(runNum,Eb,useTB,useVolatile);
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
                                if(count%10000 == 0) System.out.println(count/1000 + "k events (this is LTCC analysis on "+runstrg+") ; progress : "+progresscount+"/"+filetot);
                        }   
                        reader.close();
                }   
                System.out.println("Total : " + count + " events");
                ana.plot();
        }   
}
