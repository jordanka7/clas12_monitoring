
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

public class tof_monitor {
		boolean userTimeBased, write_volatile;
		public int runNum;
		public boolean hasRF;
		public float RFTime;
		public H2F p1a_pad_occ, p1b_pad_occ, p2_pad_occ;
		public H2F p1a_pad_XY, p1b_pad_XY, p2_pad_XY;
		public H2F[] p1a_pad_vt, p1b_pad_vt, p2_pad_vt;
		public H2F[] p1a_pad_edep, p1b_pad_edep, p2_pad_edep;
		public H2F[] p1a_pad_dt, p1b_pad_dt, p2_pad_dt;
		public H2F[][] DC_residuals_trkDoca;
		public H1F[][] DC_residuals, DC_time;
		public F1D[][] f_time_invertedS;
	public tof_monitor(int reqrunNum, boolean reqTimeBased, boolean reqwrite_volatile) {
			runNum = reqrunNum;userTimeBased=reqTimeBased;
			write_volatile = reqwrite_volatile;
			p1a_pad_occ = new H2F("p1a_pad_occ","p1a_pad_occ",25,0,25,6,0.5,6.5);
			p1a_pad_occ.setTitle("p1a occupancies");
			p1a_pad_occ.setTitleX("paddle");
			p1a_pad_occ.setTitleY("sector");
			p1b_pad_occ = new H2F("p1b_pad_occ","p1b_pad_occ",65,0,65,6,0.5,6.5);
			p1b_pad_occ.setTitle("p1b occupancies");
			p1b_pad_occ.setTitleX("paddle");
			p1b_pad_occ.setTitleY("sector");
			p2_pad_occ = new H2F("p2_pad_occ","p2_pad_occ",5,1,6,6,0.5,6.5);
			p2_pad_occ.setTitle("p2 occupancies");
			p2_pad_occ.setTitleX("paddle");
			p2_pad_occ.setTitleY("sector");
			p1a_pad_XY = new H2F("p1a_pad_XY","p1a_pad_XY",100,-500,500,100,-500,500);
			p1a_pad_XY.setTitle("p1a position");
			p1a_pad_XY.setTitleX("X (cm)");
			p1a_pad_XY.setTitleY("Y (cm)");
			p1b_pad_XY = new H2F("p1b_pad_XY","p1b_pad_XY",100,-500,500,100,-500,500);
			p1b_pad_XY.setTitle("p1b position");
			p1b_pad_XY.setTitleX("X (cm)");
			p1b_pad_XY.setTitleY("Y (cm)");
			p2_pad_XY = new H2F("p2_pad_XY","p2_pad_XY",100,-500,500,100,-500,500);
			p2_pad_XY.setTitle("p2 position");
			p2_pad_XY.setTitleX("X (cm)");
			p2_pad_XY.setTitleY("Y (cm)");

			p1a_pad_vt = new H2F[6];
			p1b_pad_vt = new H2F[6];
			p2_pad_vt = new H2F[6];
			p1a_pad_edep = new H2F[6];
			p1b_pad_edep = new H2F[6];
			p2_pad_edep = new H2F[6];
			p1a_pad_dt = new H2F[6];
			p1b_pad_dt = new H2F[6];
			p2_pad_dt = new H2F[6];
			DC_residuals_trkDoca = new H2F[6][6];
			DC_residuals = new H1F[6][6];
			DC_time = new H1F[6][6];
			f_time_invertedS = new F1D[6][6];
		for(int s=0;s<6;s++){
			p1a_pad_vt[s] = new H2F(String.format("p1a_pad_vt_S%d",s+1),String.format("p1a_pad_vt_S%d",s+1),25,0,25,100,-1.002,1.002);
			p1a_pad_vt[s].setTitle(String.format("p1a S%d time",s+1));
			p1a_pad_vt[s].setTitleX("paddle");
			p1a_pad_vt[s].setTitleY("time");
			p1b_pad_vt[s] = new H2F(String.format("p1b_pad_vt_S%d",s+1),String.format("p1b_pad_vt_S%d",s+1),65,0,65,100,-1.002,1.002);
			p1b_pad_vt[s].setTitle(String.format("p1b S%d time",s+1));
			p1b_pad_vt[s].setTitleX("paddle");
			p1b_pad_vt[s].setTitleY("time");
			p2_pad_vt[s] = new H2F(String.format("p2_pad_vt_S%d",s+1),String.format("p2_pad_vt_S%d",s+1),5,1,6,100,-1.002,1.002);
			p2_pad_vt[s].setTitle(String.format("p2 S%d time",s+1));
			p2_pad_vt[s].setTitleX("paddle");
			p2_pad_vt[s].setTitleY("time");
			p1a_pad_edep[s] = new H2F(String.format("p1a_pad_edep_S%d",s+1),String.format("p1a_pad_edep_S%d",s+1),25,0,25,100,0,50);
			p1a_pad_edep[s].setTitle(String.format("p1a S%d energy",s+1));
			p1a_pad_edep[s].setTitleX("paddle");
			p1a_pad_edep[s].setTitleY("E (MeV)");
			p1b_pad_edep[s] = new H2F(String.format("p1b_pad_edep_S%d",s+1),String.format("p1b_pad_edep_S%d",s+1),65,0,65,100,0,50);
			p1b_pad_edep[s].setTitle(String.format("p1b S%d energy",s+1));
			p1b_pad_edep[s].setTitleX("paddle");
			p1b_pad_edep[s].setTitleY("E (MeV)");
			p2_pad_edep[s] = new H2F(String.format("p2_pad_edep_S%d",s+1),String.format("p2_pad_edep_S%d",s+1),5,1,6,100,0,50);
			p2_pad_edep[s].setTitle(String.format("p2 S%d energy",s+1));
			p2_pad_edep[s].setTitleX("paddle");
			p2_pad_edep[s].setTitleY("E (MeV)");
			p1a_pad_dt[s] = new H2F(String.format("p1a_pad_dt_S%d",s+1),String.format("p1a_pad_dt_S%d",s+1),25,0,25,100,-2.004,2.004);
			p1a_pad_dt[s].setTitle(String.format("p1a S%d #delta t",s+1));
			p1a_pad_dt[s].setTitleX("paddle");
			p1a_pad_dt[s].setTitleY("time");
			p1b_pad_dt[s] = new H2F(String.format("p1b_pad_dt_S%d",s+1),String.format("p1b_pad_dt_S%d",s+1),65,0,65,100,-2.004,2.004);
			p1b_pad_dt[s].setTitle(String.format("p1b S%d #delta t",s+1));
			p1b_pad_dt[s].setTitleX("paddle");
			p1b_pad_dt[s].setTitleY("time");
			//p2_pad_dt[s] = new H2F(String.format("p2_pad_dt_S%d",s+1),String.format("p2_pad_dt_S%d",s+1),5,1,6,100,-2.004,2.004);
			p2_pad_dt[s] = new H2F(String.format("p2_pad_dt_S%d",s+1),String.format("p2_pad_dt_S%d",s+1),5,1,6,100,-12,12);
			p2_pad_dt[s].setTitle(String.format("p2 S%d #delta t",s+1));
			p2_pad_dt[s].setTitleX("paddle");
			p2_pad_dt[s].setTitleY("time");
			float[] DCcellsizeSL = {0.9f,0.9f,1.3f,1.3f,2.0f,2.0f};
			for(int sl=0;sl<6;sl++){
				DC_residuals_trkDoca[s][sl] = new H2F(String.format("DC_residuals_trkDoca_%d_%d",s+1,sl+1),String.format("DC_residuals_trkDoca_%d_%d",s+1,sl+1)
						,100,0,DCcellsizeSL[sl],100,-1,1);
				DC_residuals_trkDoca[s][sl].setTitle(String.format("DC residuals S%d SL%d",s+1,sl+1));
				DC_residuals_trkDoca[s][sl].setTitleX("DOCA (cm)");
				DC_residuals_trkDoca[s][sl].setTitleY("residual (cm)");
				DC_residuals[s][sl] = new H1F(String.format("DC_residuals_%d_%d",s+1,sl+1),String.format("DC_residuals_%d_%d",s+1,sl+1),100,-1,1);
				DC_residuals[s][sl].setTitle(String.format("DC residuals S%d SL%d",s+1,sl+1));
				DC_residuals[s][sl].setTitleX("residual (cm)");
				DC_time[s][sl] = new H1F(String.format("DC_Time_%d_%d",s+1,sl+1),String.format("DC_Time_%d_%d",s+1,sl+1),200,-100,1000);
                                DC_time[s][sl].setTitle(String.format("DC Time S%d SL%d",s+1,sl+1));
                                DC_time[s][sl].setTitleX("time (ns)");
				//DC_time[s][sl].setOptStat(1111111);
                        	DC_time[s][sl].setLineWidth(4);
				f_time_invertedS[s][sl] = new F1D(String.format("Inverted_S_%d_%d",s+1,sl+1),"[p0]/(1+exp(-[p1]*(x-[p2])))",-100,1000);
				f_time_invertedS[s][sl].setOptStat("111111");
			}
		}
	}
	/*
       {"name":"sector",       "id":3, "type":"int8",   "info":"DC sector"},
        {"name":"superlayer",   "id":4, "type":"int8",  "info":"DC superlayer (1...6)"},
        {"name":"doca",         "id":8,  "type":"float", "info":"doca of the hit calculated from TDC (in cm)"},
        {"name":"trkDoca",      "id":10,  "type":"float", "info":"track doca of the hit (in cm)"},
        {"name":"timeResidual", "id":11,  "type":"float", "info":"time residual of the hit (in cm)"},
	 */
	public void fillDC(DataBank DCB){
		for(int r=0;r<DCB.rows();r++){
			int s = DCB.getInt("sector",r)-1;
			int sl = DCB.getInt("superlayer",r)-1;
			float trkDoca = DCB.getFloat("trkDoca",r);
			float timeResidual = DCB.getFloat("timeResidual",r);
			float time = DCB.getFloat("time",r);
			float field = DCB.getFloat("B",r);
			if(s>-1&&s<6&&sl>-1&&sl<6){
				boolean otherregions = (sl<2 || sl>3);
				boolean region2 = ((sl==2||sl==3) && field<0.5);
				if (otherregions||region2) {
					DC_residuals_trkDoca[s][sl].fill(trkDoca,timeResidual);
					DC_residuals[s][sl].fill(timeResidual);
					DC_time[s][sl].fill(time);
				}
			}
			else System.out.println("sector "+(s+1)+" superlayer "+(sl+1));
		}
	}
	public void fillTOFHists(DataBank tofB, DataBank DCB){
		for(int r=0;r<tofB.rows();r++){
			int layer = tofB.getInt("layer", r);
			float thisTime = tofB.getFloat("time", r) - tofB.getFloat("pathLength", r)/29.98f - RFTime;
			thisTime = (thisTime+1.002f) % 2.004f;
			thisTime = thisTime - 1.002f;
			float thisChi2 = 999;
			float thisMom = 0;
			float thisVz = 0;
			boolean foundTrk = false;
			for(int s=0;s<DCB.rows() && !foundTrk;s++){
				//if(DCB.getInt("id",s) == tofB.getInt("trackid",r) && DCB.getInt("q",s)>0){}
				//if(DCB.getInt("id",s) == tofB.getInt("trackid",r) && DCB.getInt("q",s)<0){}
				if(DCB.getInt("id",s) == tofB.getInt("trackid",r) ){
					thisChi2 = DCB.getFloat("chi2",s);
					thisMom=0;
					thisMom += DCB.getFloat("p0_x",s)*DCB.getFloat("p0_x",s);
					thisMom += DCB.getFloat("p0_y",s)*DCB.getFloat("p0_y",s);
					thisMom += DCB.getFloat("p0_z",s)*DCB.getFloat("p0_z",s);
					thisMom = (float)Math.sqrt(thisMom);
					thisVz = DCB.getFloat("Vtx0_z",s);
					// 200 MeV
					// 400 chi2
					// 30 cm vz
					if(thisChi2 < 500 && thisMom > 0.8 && Math.abs(thisVz)<15 )foundTrk = true;
				}
			}
			if(foundTrk && tofB.getFloat("energy", r) > 1.0 ){
				//TEMPORARY TEST FOR PION VERTEX TIME
				float flightTime = tofB.getFloat("pathLength", r)/(float)( 29.98f * thisMom/Math.sqrt(thisMom*thisMom + 0.13957f*0.13957f) );
				float thisPionTime = tofB.getFloat("time", r) - flightTime - RFTime;
				thisPionTime = (thisPionTime+1.002f) % 2.004f;
				thisPionTime = thisPionTime - 1.002f;
				//System.out.println("compare vextex times : "+thisTime+" =? "+thisPionTime+" , diff = "+(thisTime-thisPionTime));
				//thisTime = thisPionTime;
				thisVz = 0;
				if(layer==1){//panel 1A
					p1a_pad_occ.fill( tofB.getInt("component", r) , tofB.getInt("sector", r) );
					p1a_pad_XY.fill( tofB.getFloat("x",r) , tofB.getFloat("y",r) );
					p1a_pad_vt[tofB.getInt("sector", r)-1].fill( tofB.getInt("component", r) , thisTime-thisVz/29.98f );
					p1a_pad_edep[tofB.getInt("sector", r)-1].fill( tofB.getInt("component", r) , tofB.getFloat("energy", r) );
				}
				if(layer==2){//panel 1B
					p1b_pad_occ.fill( tofB.getInt("component", r) , tofB.getInt("sector", r) );
					p1b_pad_XY.fill( tofB.getFloat("x",r) , tofB.getFloat("y",r) );
					p1b_pad_vt[tofB.getInt("sector", r)-1].fill( tofB.getInt("component", r) , thisTime-thisVz/29.98f );
					p1b_pad_edep[tofB.getInt("sector", r)-1].fill( tofB.getInt("component", r) , tofB.getFloat("energy", r) );
				}
				if(layer==3){//panel 2
					p2_pad_occ.fill( tofB.getInt("component", r) , tofB.getInt("sector", r) );
					p2_pad_XY.fill( tofB.getFloat("x",r) , tofB.getFloat("y",r) );
					p2_pad_vt[tofB.getInt("sector", r)-1].fill( tofB.getInt("component", r) , thisTime-thisVz/29.98f );
					p2_pad_edep[tofB.getInt("sector", r)-1].fill( tofB.getInt("component", r) , tofB.getFloat("energy", r) );
				}
			}
		}
		for(int h1=0;h1<tofB.rows();h1++){
			float thisChi2 = 999;
			float thisMom = 0;
			boolean foundTrk1 = false;int iTrk1=-1;
			boolean foundTrk2 = false;int iTrk2=-1;
			if(tofB.getInt("layer", h1)==1&&tofB.getFloat("energy", h1)>1.0)for(int t=0;t<DCB.rows() && !foundTrk1;t++){
				if(DCB.getInt("id",t) == tofB.getInt("trackid",h1) ){
					thisChi2 = DCB.getFloat("chi2",t);
					thisMom=0;
					thisMom += DCB.getFloat("p0_x",t)*DCB.getFloat("p0_x",t);
					thisMom += DCB.getFloat("p0_y",t)*DCB.getFloat("p0_y",t);
					thisMom += DCB.getFloat("p0_z",t)*DCB.getFloat("p0_z",t);
					thisMom = (float)Math.sqrt(thisMom);
					if(thisChi2 < 5000 && thisMom > 0.8 && Math.abs(DCB.getFloat("Vtx0_z",t))<15 ){foundTrk1 = true;iTrk1=t;}
				}
			}
			if(foundTrk1)for(int h2=0;h2<tofB.rows();h2++){
				if(h1!=h2 && tofB.getInt("layer", h2)!=1 && tofB.getInt("sector", h2)!=tofB.getInt("sector", h1) && tofB.getFloat("energy", h2)>1.0 )
				  for(int t=0;t<DCB.rows() && !foundTrk2;t++){
					if(t!=iTrk1 && DCB.getInt("id",t) == tofB.getInt("trackid",h2) ){
						thisChi2 = DCB.getFloat("chi2",t);
						thisMom=0;
						thisMom += DCB.getFloat("p0_x",t)*DCB.getFloat("p0_x",t);
						thisMom += DCB.getFloat("p0_y",t)*DCB.getFloat("p0_y",t);
						thisMom += DCB.getFloat("p0_z",t)*DCB.getFloat("p0_z",t);
						thisMom = (float)Math.sqrt(thisMom);
						if(thisChi2 < 5000 && thisMom > 0.8 && Math.abs(DCB.getFloat("Vtx0_z",t))<15 ){foundTrk2 = true;iTrk2=t;}
					}
				}
				if(foundTrk2){
					float thisTime = tofB.getFloat("time", h1) - tofB.getFloat("time", h2) - (tofB.getFloat("pathLength",h1)-tofB.getFloat("pathLength",h2))/29.98f;;
					if(thisTime!=0){
						p1a_pad_dt[tofB.getInt("sector",h1)-1].fill( tofB.getInt("component", h1) , thisTime );
						if(tofB.getInt("layer", h2)==2)p1b_pad_dt[tofB.getInt("sector",h2)-1].fill( tofB.getInt("component", h2) , thisTime );
						if(tofB.getInt("layer", h2)==3)p2_pad_dt[tofB.getInt("sector",h2)-1].fill( tofB.getInt("component", h2) , thisTime );
					}
				}
			}
		}
	}
	public void fillRFTime(DataBank RFB){
		for(int r=0;r<RFB.rows() && !hasRF;r++){
			if(RFB.getInt("id",r)==1){
				hasRF=true;
				RFTime = RFB.getFloat("time",r);
			}
		}
	}
        public void processEvent(DataEvent event) {
		hasRF = false;
		DataBank trackDetBank = null, hitBank = null;
		if(userTimeBased){
			if(event.hasBank("TimeBasedTrkg::TBTracks"))trackDetBank = event.getBank("TimeBasedTrkg::TBTracks");
			if(event.hasBank("TimeBasedTrkg::TBHits")){hitBank = event.getBank("TimeBasedTrkg::TBHits");}
		}
		if(!userTimeBased){
			if(event.hasBank("HitBasedTrkg::HBTracks"))trackDetBank = event.getBank("HitBasedTrkg::HBTracks");
			if(event.hasBank("HitBasedTrkg::HBHits"))hitBank = event.getBank("HitBasedTrkg::HBHits");
		}
		if(event.hasBank("RUN::rf"))fillRFTime(event.getBank("RUN::rf"));
		if(!hasRF)return;
		if(event.hasBank("FTOF::hits") && trackDetBank!=null)fillTOFHists(event.getBank("FTOF::hits") , trackDetBank);
		if(userTimeBased && hitBank!=null)fillDC(hitBank);

	}

	public void initInvertedSFitPar(int slayer, F1D function) {
		double min = 100.;
		double max = 220.;
		if (slayer == 1) {
			min = 100.; max = 240.;
			function.setParameter(1,-0.038); function.setParLimits(1,-0.01,-0.06);
			function.setParameter(2,118.); function.setParLimits(2,100.,150.);
		}
		if (slayer == 2) {
			min = 120.; max = 240.;
			function.setParameter(1,-0.040); function.setParLimits(1,-0.01,-0.06);
			function.setParameter(2,136.); function.setParLimits(2,100.,200.);
		}
		if (slayer == 3) {
			min = 200.; max = 450.;
			function.setParameter(1,-0.030);function.setParLimits(1,-0.01,-0.05);
			function.setParameter(2,320.); function.setParLimits(2,200.,500.);
		}
		if (slayer == 4) {
			min = 200.; max = 500.;
			function.setParameter(1,-0.023); function.setParLimits(1,-0.01,-0.05);
			function.setParameter(2,350.); function.setParLimits(2,200.,500.);
		}
		if (slayer == 5) {
			min = 400.; max = 700.;
			function.setParameter(1,-0.024);function.setParLimits(1,-0.01,-0.05);
			function.setParameter(2,623.); function.setParLimits(2,500.,700.);
		}
		if (slayer == 6) {
			min = 480.; max = 700.;
			function.setParameter(1,-0.034); function.setParLimits(1,-0.01,-0.05);
			function.setParameter(2,683.); function.setParLimits(2,500.,750.);
		}
		function.setRange(min,max);
		function.setLineColor(2);
		function.setLineWidth(4);
	}

	public void analyze() {
		for(int sl=1;sl<=6;sl++) {
			for (int l=0;l<6;l++) {
				initInvertedSFitPar(sl,f_time_invertedS[l][sl-1]);
				DataFitter.fit(f_time_invertedS[l][sl-1],DC_time[l][sl-1],"QL");
				//System.out.println("Fit DC Time for sector "+l+ " superlayer " +sl+ "complete");
				//System.out.println(" ");System.out.println(" ");
				DC_time[l][sl-1].setFunction(null);
			}
		}
	}

        public void plot() {
                EmbeddedCanvas can_TOF_occ = new EmbeddedCanvas();
                can_TOF_occ.setSize(3000,5000);
                can_TOF_occ.divide(6,9);
                can_TOF_occ.setAxisTitleSize(18);
                can_TOF_occ.setAxisFontSize(18);
                can_TOF_occ.setTitleSize(18);
                for(int s=0;s<6;s++){
			can_TOF_occ.cd(s);can_TOF_occ.draw(p1a_pad_vt[s]);
			//can_TOF_occ.getPad(6+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(6+s);can_TOF_occ.draw(p1b_pad_vt[s]);
			//can_TOF_occ.getPad(12+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(12+s);can_TOF_occ.draw(p2_pad_vt[s]);
			//can_TOF_occ.getPad(12+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(18+s);can_TOF_occ.draw(p1a_pad_edep[s]);
			//can_TOF_occ.getPad(24+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(24+s);can_TOF_occ.draw(p1b_pad_edep[s]);
			//can_TOF_occ.getPad(30+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(30+s);can_TOF_occ.draw(p2_pad_edep[s]);
			//can_TOF_occ.getPad(36+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(36+s);can_TOF_occ.draw(p1a_pad_dt[s]);
			//can_TOF_occ.getPad(42+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(42+s);can_TOF_occ.draw(p1b_pad_dt[s]);
			//can_TOF_occ.getPad(48+s).getAxisZ().setLog(true);
			can_TOF_occ.cd(48+s);can_TOF_occ.draw(p2_pad_dt[s]);
			//can_TOF_occ.getPad(54+s).getAxisZ().setLog(true);
		}
		if(runNum>0){
			if(!write_volatile)can_TOF_occ.save(String.format("plots"+runNum+"/TOF_cal.png"));
			if(write_volatile)can_TOF_occ.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/TOF_cal.png"));
			System.out.println(String.format("saved plots"+runNum+"/TOF_cal.png"));
		}
		else{
			can_TOF_occ.save(String.format("plots/TOF_cal.png"));
			System.out.println(String.format("saved plots/TOF_cal.png"));
		}

		EmbeddedCanvas can_DC_resd_trkDoca  = new EmbeddedCanvas();
		can_DC_resd_trkDoca.setSize(3000,3000);
		can_DC_resd_trkDoca.divide(6,6);
		can_DC_resd_trkDoca.setAxisTitleSize(18);
		can_DC_resd_trkDoca.setAxisFontSize(18);
		can_DC_resd_trkDoca.setTitleSize(18);
		for(int s=0;s<6;s++)for(int sl=0;sl<6;sl++){
			can_DC_resd_trkDoca.cd(sl + 6*s);
			can_DC_resd_trkDoca.getPad(sl + 6*s).getAxisZ().setLog(true);
			can_DC_resd_trkDoca.draw(DC_residuals_trkDoca[s][sl]);
		}
		if(runNum>0){
			if(!write_volatile)can_DC_resd_trkDoca.save(String.format("plots"+runNum+"/DC_resd_trkDoca.png"));
			if(write_volatile)can_DC_resd_trkDoca.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/DC_resd_trkDoca.png"));
			System.out.println(String.format("saved plots"+runNum+"/DC_resd_trkDoca.png"));
		}
		else{
			can_DC_resd_trkDoca.save(String.format("plots/DC_resd_trkDoca.png"));
			System.out.println(String.format("saved plots/DC_resd_trkDoca.png"));
		}

		EmbeddedCanvas can_DC_resd  = new EmbeddedCanvas();
		can_DC_resd.setSize(3000,3000);
		can_DC_resd.divide(6,6);
		can_DC_resd.setAxisTitleSize(18);
		can_DC_resd.setAxisFontSize(18);
		can_DC_resd.setTitleSize(18);
		for(int s=0;s<6;s++)for(int sl=0;sl<6;sl++){
			can_DC_resd.cd(sl + 6*s);
			can_DC_resd.draw(DC_residuals[s][sl]);
		}
		if(runNum>0){
			if(!write_volatile)can_DC_resd.save(String.format("plots"+runNum+"/DC_resd.png"));
			if(write_volatile)can_DC_resd.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/DC_resd.png"));
			System.out.println(String.format("saved plots"+runNum+"/DC_resd.png"));
		}
		else{
			can_DC_resd.save(String.format("plots/DC_resd.png"));
			System.out.println(String.format("saved plots/DC_resd.png"));
		}

		EmbeddedCanvas can_DC_time  = new EmbeddedCanvas();
                can_DC_time.setSize(3000,3000);
                can_DC_time.divide(6,6);
                can_DC_time.setAxisTitleSize(18);
                can_DC_time.setAxisFontSize(18);
                can_DC_time.setTitleSize(18);
                for(int s=0;s<6;s++)for(int sl=0;sl<6;sl++){
                        can_DC_time.cd(sl + 6*s);
                        can_DC_time.draw(DC_time[s][sl]);can_DC_time.draw(f_time_invertedS[s][sl],"same");
                }
                if(runNum>0){
                        if(!write_volatile)can_DC_time.save(String.format("plots"+runNum+"/DC_time.png"));
                        if(write_volatile)can_DC_time.save(String.format("/volatile/clas12/rga/spring18/plots"+runNum+"/DC_time.png"));
                        System.out.println(String.format("saved plots"+runNum+"/DC_time.png"));
                }
                else{
                        can_DC_time.save(String.format("plots/DC_time.png"));
                        System.out.println(String.format("saved plots/DC_time.png"));
                }
	}

        public void write() {
                TDirectory dirout = new TDirectory();
                dirout.mkdir("/tof/");
                dirout.cd("/tof/");
                for(int s=0;s<6;s++){
                        dirout.addDataSet(p1a_pad_vt[s],p1b_pad_vt[s],p2_pad_vt[s],p1a_pad_dt[s],p1b_pad_dt[s],p2_pad_dt[s]);
		}
		dirout.mkdir("/dc/");
		dirout.cd("/dc/");
		for(int s=0;s<6;s++)for(int sl=0;sl<6;sl++){
			dirout.addDataSet(DC_residuals_trkDoca[s][sl],DC_time[s][sl]);
		}
                if(write_volatile)if(runNum>0)dirout.writeFile("/volatile/clas12/rga/spring18/plots"+runNum+"/out_TOF_"+runNum+".hipo");

		if(!write_volatile){
			if(runNum>0)dirout.writeFile("plots"+runNum+"/out_TOF_"+runNum+".hipo");
			else dirout.writeFile("plots/out_TOF.hipo");
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
		if(args.length>0)runNum = Integer.parseInt(args[0]);
		if(args.length>1)filelist = args[1];
		if(args.length>2)if(Integer.parseInt(args[2])==0)useTB=false;
		tof_monitor ana = new tof_monitor(runNum,useTB,useVolatile);
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
		int maxevents = 50000000;
		if(args.length>2)maxevents=Integer.parseInt(args[2]);
		int progresscount=0;int filetot = toProcessFileNames.size();
		for (String runstrg : toProcessFileNames) if(count<maxevents){
			progresscount++;
			System.out.println(String.format(">>>>>>>>>>>>>>>> tof_monitor %s",runstrg));
			File varTmpDir = new File(runstrg);
			if(!varTmpDir.exists()){System.out.println("FILE DOES NOT EXIST");continue;}
			System.out.println("READING NOW "+runstrg);
			HipoDataSource reader = new HipoDataSource();
			reader.open(runstrg);
			int filecount = 0;
			while(reader.hasEvent() && count<maxevents) {
				DataEvent event = reader.getNextEvent();
				ana.processEvent(event);
				filecount++;count++;
				if(count%10000 == 0) System.out.println(count/1000 + "k events (this is tof_monitor on "+runstrg+") progress : " + progresscount+"/"+filetot);
			}
			reader.close();
		}
		System.out.println("Total : " + count + " events");
		ana.analyze();
		ana.plot();
		ana.write();
        }
}
