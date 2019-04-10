
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

public class central {
	boolean userTimeBased, write_volatile;
	public int runNum;
	public boolean BackToBack;
	public float STT, MinCTOF,MaxCTOF, minSTT, maxSTT;
	public float[] CTOF_shft;
	public H2F H_CTOF_pos, H_CTOF_edep_phi, H_CTOF_edep_z, H_CTOF_path_mom;
	public H2F H_vz_DC_CVT, H_phi_DC_CVT, H_CVT_CTOF_phi, H_CVT_CTOF_z, H_CVT_t_STT, H_CVT_t_pad;
	public H1F[] H_CVT_t;
	public central(int reqrunNum, boolean reqTimeBased, boolean reqwrite_volatile) {
		runNum = reqrunNum;userTimeBased=reqTimeBased;
		write_volatile = reqwrite_volatile;
		//CTOF_shft = new float[]{0.00f , -108.81f , -105.97f , -103.48f , -109.72f , -107.39f , -104.89f , -109.08f , -107.77f , -104.47f , 
		//			-109.55f , -107.23f , -105.16f , -108.38f , 0.00f , -104.24f , -108.10f , -106.98f , -104.59f , -108.12f , 
		//			-106.89f , -104.63f , -108.99f , -106.41f , -105.22f , -108.29f , -106.25f , -103.80f , -108.60f , -106.13f , 
		//			-103.67f , -107.81f , -106.08f , -103.66f , -108.02f , -106.51f , -105.12f , -108.38f , -107.36f , -104.17f , 
		//			-108.31f , -105.56f , -102.91f , -107.76f , -106.62f , -103.89f , -108.73f , -106.88f , -104.80f , -106.36f};//2052

		//CTOF_shft = new float[]{0.00f , -557.95f , -552.60f , -560.09f , -556.33f , -557.89f , -559.68f , -559.52f , -563.12f , -559.26f , 
		//			-555.64f , -555.37f , -566.83f , -558.37f , 0.00f , -557.06f , -572.78f , -205.33f , -568.07f , -200.87f , 
		//			-550.44f , -222.76f , -541.37f , -252.70f , -545.07f , -553.92f , -238.36f , -554.90f , -218.27f , -539.17f , 
		//			-224.08f , -543.47f , -222.41f , -553.31f , -542.48f , -542.28f , -557.34f , -541.96f , -552.40f , -551.26f , 
		//			-542.06f , -239.25f , -539.38f , -229.84f , -563.36f , -266.99f , -549.68f , -203.33f , -555.38f , -420.01f};//2195
		//CTOF_shft = new float[]{0.00f , -299.73f , -297.20f , -294.74f , -300.75f , -298.03f , -295.33f , -299.77f , -298.33f , -294.47f , 
		//			-300.32f , -297.32f , -295.15f , -299.03f , 0.00f , -294.54f , -298.85f , -298.03f , -295.60f , -299.52f , 
		//			-297.84f , -294.92f , -299.51f , -296.78f , -295.84f , -299.06f , -296.87f , -294.73f , -299.62f , -296.65f , 
		//			-294.08f , -298.50f , -296.86f , -294.34f , -299.30f , -297.71f , -295.99f , -299.43f , -298.16f , -294.76f , 
		//			-299.15f , -296.72f , -294.17f , -299.01f , -297.59f , -294.53f , -299.59f , -297.51f , -295.47f , -297.36f};//2383
		//CTOF_shft = new float[]{0.00f , -301.09f , -298.33f , -295.57f , -301.58f , -298.70f , -295.93f , -300.22f , -298.70f , -295.47f , 
		//			-300.80f , -297.82f , -295.74f , -299.63f , 0.00f , -294.88f , -299.22f , -298.66f , -295.95f , -299.18f , 
		//			-297.50f , -294.83f , -299.31f , -296.56f , -295.67f , -299.12f , -296.83f , -294.34f , -299.45f , -296.44f , 
		//			-293.88f , -298.34f , -296.83f , -295.05f , -299.40f , -297.57f , -295.61f , -299.46f , -298.14f , -294.77f , 
		//			-299.22f , -297.05f , -294.38f , -299.45f , -298.18f , -295.27f , -300.25f , -298.39f , -296.41f , -297.84f};//2478 
                CTOF_shft = new float[]{0.00f , -300.51f , -297.72f , -295.05f , -301.03f , -298.19f , -295.41f , -299.74f , -298.30f , -295.21f ,
                                        -300.22f , -297.35f , -295.28f , -299.14f , 0.00f , -294.51f , -298.80f , -298.30f , -295.49f , -298.81f ,
                                        -297.23f , -294.44f , -298.96f , -296.27f , -295.30f , -298.74f , -296.44f , -293.98f , -299.07f , -296.08f ,
                                        -293.52f , -297.96f , -296.42f , -294.40f , -298.90f , -297.14f , -295.28f , -299.03f , -297.71f , -294.35f ,
                                        -298.81f , -296.75f , -293.77f , -298.96f , -297.65f , -294.67f , -299.73f , -297.76f , -295.87f , -297.10f};//2476

		H_CTOF_pos = new H2F("H_CTOF_pos","H_CTOF_pos",50,-180,180,100,-10,5);
		//H_CTOF_pos = new H2F("H_CTOF_pos","H_CTOF_pos",50,-180,180,100,-20,10);
		H_CTOF_pos.setTitle("CTOF hits z vs #phi");
		H_CTOF_pos.setTitleX("#phi (^o)");
		H_CTOF_pos.setTitleY("z (cm)");
		H_CTOF_edep_phi = new H2F("H_CTOF_edep_phi","H_CTOF_edep_phi",50,-180,180,100,0,150);
		H_CTOF_edep_phi.setTitle("CTOF Edep vs #phi");
		H_CTOF_edep_phi.setTitleX("#phi (^o)");
		H_CTOF_edep_phi.setTitleY("E (MeV)");
		H_CTOF_edep_z = new H2F("H_CTOF_edep_z","H_CTOF_edep_z",100,-10,5,100,0,150);
		//H_CTOF_edep_z = new H2F("H_CTOF_edep_z","H_CTOF_edep_z",100,-20,10,100,0,150);
		H_CTOF_edep_z.setTitle("CTOF Edep vs z");
		H_CTOF_edep_z.setTitleX("z (cm)");
		H_CTOF_edep_z.setTitleY("E (MeV)");
		H_vz_DC_CVT = new H2F("H_vz_DC_CVT","H_vz_DC_CVT",100,-20,20,100,-20,20);
		H_vz_DC_CVT.setTitle("CVT vz vs DC vz");
		H_vz_DC_CVT.setTitleX("DC vz (cm)");
		H_vz_DC_CVT.setTitleY("CVT vz (cm)");
		H_phi_DC_CVT = new H2F("H_phi_DC_CVT","H_phi_DC_CVT",100,-180,180,100,-180,180);
		H_phi_DC_CVT.setTitle("CVZ #phi vs DC #phi");
		H_phi_DC_CVT.setTitleX("DC #phi (^o)");
		H_phi_DC_CVT.setTitleY("CVT #phi (^o)");
		H_CTOF_path_mom = new H2F("H_CTOF_path_mom","H_CTOF_path_mom",100,0,2,100,20,60);
		H_CTOF_path_mom.setTitle("CTOT path vs mom");
		H_CTOF_path_mom.setTitleX("p (GeV)");
		H_CTOF_path_mom.setTitleY("CTOF path (cm)");
		H_CVT_CTOF_phi = new H2F("H_CVT_CTOF_phi","H_CVT_CTOF_phi",100,-180,180,50,-180,180);
		H_CVT_CTOF_phi.setTitle("CVT CTOF #Delta #phi");
		H_CVT_CTOF_phi.setTitleX("CVT #phi (^o)");
		H_CVT_CTOF_phi.setTitleY("CTOF #phi (^o)");
		H_CVT_CTOF_z = new H2F("H_CVT_CTOF_z","H_CVT_CTOF_z",100,-10,5,100,-10,5);
		H_CVT_CTOF_z.setTitle("CVT CTOF #Delta z");
		H_CVT_CTOF_z.setTitleX("CVT z (cm)");
		H_CVT_CTOF_z.setTitleY("CTOF z (cm)");
		minSTT = 100;maxSTT=300;
		if(runNum>0&&runNum<3210){
			//minSTT = 150;
			//maxSTT = 250;
			minSTT = 540;maxSTT=590;
		}
		//MinCTOF=-5;MaxCTOF=5;
		MinCTOF=320;MaxCTOF=400;
		if(runNum>0&&runNum<3210){
			MinCTOF=-5;MaxCTOF=5;
		}
		MinCTOF=-5;MaxCTOF=5;
		H_CVT_t_STT = new H2F("H_CVT_t_STT","H_CVT_t_STT",100,minSTT,maxSTT,100,minSTT+(MinCTOF+MaxCTOF)/2,maxSTT+(MinCTOF+MaxCTOF)/2);
		//H_CVT_t_STT = new H2F("H_CVT_t_STT","H_CVT_t_STT",100,150,250,100,520,620);
		H_CVT_t_STT.setTitle("CTOF vertex t vs STT");
		H_CVT_t_STT.setTitleX("STT (ns)");
		H_CVT_t_STT.setTitleY("CTOF t (ns)");
                try {
                        Thread.sleep(5000);// in ms
                }catch (Exception e) {
                        System.out.println(e);
                }
		System.out.println("CTOF range "+MinCTOF+" to "+MaxCTOF);
		//MinCTOF=-125;MaxCTOF=-75;
		//MinCTOF=-5;MaxCTOF=5;
		H_CVT_t_pad = new H2F("H_CVT_t_pad","H_CVT_t_pad",50,0.5,50.5,250,MinCTOF,MaxCTOF);
		//H_CVT_t_pad = new H2F("H_CVT_t_pad","H_CVT_t_pad",50,0.5,50.5,100,0,10);
		//H_CVT_t_pad = new H2F("H_CVT_t_pad","H_CVT_t_pad",50,0.5,50.5,100,-115,-95);//2052
		//H_CVT_t_pad = new H2F("H_CVT_t_pad","H_CVT_t_pad",50,0.5,50.5,100,-5,5);
		H_CVT_t_pad.setTitle("CVT t vs pad");
		H_CVT_t_pad.setTitleX("pad");
		H_CVT_t_pad.setTitleY("t (ns)");
		H_CVT_t = new H1F[50];
		for(int p=0;p<49;p++){
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,CTOF_shft[p]-12.5,CTOF_shft[p]+12.5);
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,CTOF_shft[p]-7.5,CTOF_shft[p]+7.5);
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,CTOF_shft[p]-2.5,CTOF_shft[p]+2.5);
			H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,MinCTOF,MaxCTOF);
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),100,0,10);
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),100,-2.5,2.5);
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),100,-115,-95);
			H_CVT_t[p].setTitle(String.format("pad %d time",p+1));
			H_CVT_t[p].setTitleX("t (ns)");
		}
		for(int p=49;p<50;p++){
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,-115,-95);
			//H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,-5,5);
			H_CVT_t[p] = new H1F(String.format("H_CVT_t_p%d",p+1),String.format("H_CVT_t_p%d",p+1),250,MinCTOF,MaxCTOF);
			H_CVT_t[p].setTitle(String.format("pad %d time",p+1));
			H_CVT_t[p].setTitleX("t (ns)");
		}
	}
	public double Vangle(Vector3 v1, Vector3 v2){ 
		double res = 0; 
		double l1 = v1.mag();
		double l2 = v2.mag();
		double prod = v1.dot(v2);
		if( l1 * l2 !=0 && Math.abs(prod)<l1*l2 )res = Math.toDegrees( Math.acos(prod/(l1*l2) ) ); 
		return res; 
	}
	public void FillTracks(DataBank DCbank, DataBank CVTbank){
		for(int iDC=0;iDC<DCbank.rows();iDC++){
			for(int iCVT=0;iCVT<CVTbank.rows() && !BackToBack;iCVT++){
				float DelPhi = (float)Math.toDegrees( Math.atan2(DCbank.getFloat("p0_y",iDC),DCbank.getFloat("p0_x",iDC)) - CVTbank.getFloat("phi0",iCVT));
				//DelPhi += 180f - 10f;
				while(DelPhi>180)DelPhi-=360;
				while(DelPhi<-180)DelPhi+=360;
				if(Math.abs(DelPhi)<180){
					H_vz_DC_CVT.fill(DCbank.getFloat("Vtx0_z",iDC),CVTbank.getFloat("z0",iCVT));
					//H_vz_DC_CVT.fill(DCbank.getFloat("Vtx0_z",iDC),0.1f*CVTbank.getFloat("z0",iCVT));
					H_phi_DC_CVT.fill(Math.toDegrees(Math.atan2(DCbank.getFloat("p0_y",iDC),DCbank.getFloat("p0_x",iDC))),Math.toDegrees(CVTbank.getFloat("phi0",iCVT)));
					BackToBack = true;
				}
			}
		}
	}
	public void FillCVTCTOF(DataBank CVTbank, DataBank CTOFbank){
		for(int iCTOF=0;iCTOF<CTOFbank.rows();iCTOF++){
			float e = CTOFbank.getFloat("energy",iCTOF);
			int trackid = CTOFbank.getInt("trkID",iCTOF);
			if(!Float.isNaN(e) && trackid>-1 && e>1.5){
				boolean matched = false;
				for(int iCVT=0;iCVT<CVTbank.rows() && !matched;iCVT++){
					float mom = CVTbank.getFloat("p",iCVT);
					float cx = CVTbank.getFloat("c_x",iCVT)*0.1f;
					float cy = CVTbank.getFloat("c_y",iCVT)*0.1f;
					float cz = CVTbank.getFloat("c_z",iCVT)*0.1f;
					float cphi = (float)Math.toDegrees(Math.atan2(cy,cx));
					int pad = CTOFbank.getInt("component",iCTOF);
					float x = CTOFbank.getFloat("x",iCTOF)*0.1f;
					float y = CTOFbank.getFloat("y",iCTOF)*0.1f;
					float z = CTOFbank.getFloat("z",iCTOF)*0.1f;
					float t = CTOFbank.getFloat("time",iCTOF);
					float p = CTOFbank.getFloat("pathLength",iCTOF);
					float phi = (float)Math.toDegrees(Math.atan2(y,x));
					float beta =  mom/(float)Math.sqrt(mom*mom+0.93827f*0.93827f);
					//float DelPhi = phi-cphi+190;
					float DelPhi = phi-cphi;
					while(DelPhi>180)DelPhi-=360;
					while(DelPhi<-180)DelPhi+=360;
					if( Math.abs(DelPhi)<30 && z<4.7){}
					if( Math.abs(DelPhi)<180){
						H_CVT_CTOF_phi.fill(cphi,phi);
						H_CVT_CTOF_z.fill(cz,z);
						H_CTOF_pos.fill(phi,z);
						H_CTOF_edep_phi.fill(phi,e);
						H_CTOF_edep_z.fill(z,e);
						H_CTOF_path_mom.fill(mom,p);
						float CTOFTime = t - p/29.92f/beta;
						//float CTOFTime = t - p/29.92f/beta - CTOF_shft[pad];
						H_CVT_t_STT.fill(STT,t - p/29.92f/beta);
						H_CVT_t_pad.fill(pad,CTOFTime-STT);
						H_CVT_t[pad].fill(CTOFTime-STT);
						H_CVT_t[49].fill(CTOFTime-STT);
						matched = true;
					}
				}
			}
		}
	}
	public void processEvent(DataEvent event) {
		BackToBack = false;
		DataBank eventBank = null, trackDetBank = null;
                if(userTimeBased){
                        if(event.hasBank("REC::Event"))eventBank = event.getBank("REC::Event");
                        if(event.hasBank("TimeBasedTrkg::TBTracks"))trackDetBank = event.getBank("TimeBasedTrkg::TBTracks");
                }
                if(!userTimeBased){
                        if(event.hasBank("RECHB::Event"))eventBank = event.getBank("RECHB::Event");
                        if(event.hasBank("HitBasedTrkg::HBTracks"))trackDetBank = event.getBank("HitBasedTrkg::HBTracks");
                }

		if(eventBank!=null)STT = eventBank.getFloat("STTime",0);
		else return;
		if(trackDetBank != null && event.hasBank("CVTRec::Tracks"))FillTracks(trackDetBank,event.getBank("CVTRec::Tracks"));
		if(BackToBack && event.hasBank("CVTRec::Tracks") && event.hasBank("CTOF::hits"))FillCVTCTOF(event.getBank("CVTRec::Tracks"),event.getBank("CTOF::hits"));
	}
	public void plot() {
		EmbeddedCanvas can_central  = new EmbeddedCanvas();
                can_central.setSize(2000,3000);
                can_central.divide(4,6);
                can_central.setAxisTitleSize(18);
                can_central.setAxisFontSize(18);
                can_central.setTitleSize(18);
		can_central.cd(0);can_central.draw(H_CTOF_pos);
		can_central.cd(1);can_central.draw(H_CTOF_edep_phi);
		can_central.cd(2);can_central.draw(H_CTOF_edep_z);
		can_central.cd(3);can_central.draw(H_CTOF_path_mom);
		can_central.cd(4);can_central.draw(H_CVT_CTOF_phi);
		can_central.cd(5);can_central.draw(H_CVT_CTOF_z);
		can_central.cd(6);can_central.draw(H_vz_DC_CVT);
		can_central.cd(7);can_central.draw(H_phi_DC_CVT);
		can_central.cd(8);can_central.draw(H_CVT_t_STT);
		can_central.cd(9);can_central.draw(H_CVT_t_pad);
		can_central.cd(10);can_central.draw(H_CVT_t[1]);for(int p=1;p<49;p++)can_central.draw(H_CVT_t[p],"same");
		can_central.getPad(10).getAxisX().setRange(MinCTOF,MaxCTOF);
		can_central.cd(11);can_central.draw(H_CVT_t[49]);
		for(int p=0;p<12;p++){
			can_central.cd(12+p);can_central.draw(H_CVT_t[16+p]);
		}
		if(runNum>0){
			if(!write_volatile)can_central.save(String.format("plots"+runNum+"/central.png"));
			if(write_volatile)can_central.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/central.png"));
			System.out.println(String.format("saved plots"+runNum+"/central.png"));
		}
		else{
			can_central.save(String.format("plots/central.png"));
			System.out.println(String.format("saved plots/central.png"));
		}
		for(int p=0;p<10;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=10;p<20;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=20;p<30;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=30;p<40;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
		for(int p=40;p<50;p++)System.out.print(String.format("%1.2ff , ",H_CVT_t[p].getMean()));
		System.out.print("\n");
	}
        public void write(){
                TDirectory dirout = new TDirectory();
                dirout.mkdir("/ctof/");
                dirout.cd("/ctof/");
                dirout.addDataSet(H_CVT_t_pad,H_CTOF_edep_phi);
                for(int p=0;p<50;p++)dirout.addDataSet(H_CVT_t[p]);
                
		if(write_volatile)if(runNum>0)dirout.writeFile("/volatile/clas12/rga/spring18/plots"+runNum+"/out_CTOF_"+runNum+".hipo");
                
		if(!write_volatile){
			if(runNum>0)dirout.writeFile("plots"+runNum+"/out_CTOF_"+runNum+".hipo");
			else dirout.writeFile("plots/out_CTOF.hipo");
		}
        }   

////////////////////////////////////////////////
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
		if(args.length>2)if(Integer.parseInt(args[2])==0)useTB=false;
		central ana = new central(runNum,useTB,useVolatile);
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
		//int maxevents = 50000000;
		int maxevents = 500000;
		if(args.length>2)maxevents=Integer.parseInt(args[2]);
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
				if(count%10000 == 0) System.out.println(count/1000 + "k events (this is central analysis on "+runstrg+") ; progress : "+progresscount+"/"+filetot);
			}
			reader.close();
		}
		System.out.println("Total : " + count + " events");
		ana.plot();
		ana.write();
        }
}
