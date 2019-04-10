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
import org.clas.viewer.AnalysisMonitor;
import org.jlab.clas.pdg.PhysicsConstants;
import org.jlab.clas.physics.Particle;
import org.jlab.groot.group.DataGroup;

public class FT {
	boolean userTimeBased, write_volatile;
	boolean[] trigger_bits;
	public int runNum, trigger;
	public boolean hasRF;
	public double startTime, rfTime;

	public double rfPeriod;
	
	public H1F summary;
	//Hodoscope
               public H1F[] hi_hodo_eall, hi_hodo_ematch, hi_hodo_tmatch;
               public F1D[] f_charge_landau;
               public H2F[] hi_hodo_ematch_2D, hi_hodo_tmatch_2D;

	//Calorimeter
	public H1F hi_cal_nclusters, hi_cal_clsize, hi_cal_clsize_ch, hi_cal_e_all, hi_cal_e_ch, hi_cal_e_neu, hi_cal_theta_ch, hi_cal_phi_ch, hi_cal_time_ch, hi_cal_time_cut_ch, hi_cal_time_neu, hi_cal_time_cut_neu; 
	public H2F hi_cal_clsize_en, hi_cal_time_e_ch, hi_cal_time_theta_ch, hi_cal_time_e_neu, hi_cal_time_theta_neu;
	public F1D ftime_ch, ftime_neu;
		
	//pi0
	public H1F hpi0sum;
	public F1D fpi0;
	public H2F hmassangle;
		
		
	public FT(int reqrunNum, boolean reqTimeBased, boolean reqwrite_volatile) {
		runNum = reqrunNum;
		userTimeBased=reqTimeBased;
		write_volatile = reqwrite_volatile;

		trigger_bits = new boolean[32];
		startTime=-1000;
		rfTime=-1000;
		trigger = 0;
		
		rfPeriod = 4.008;

		H1F summary = new H1F("summary","summary",6,1,7);
       		summary.setTitleX("sector");
       		summary.setTitleY("DC hits");
       		summary.setFillColor(33);

		//Hodoscope Histograms
		hi_hodo_eall = new H1F[2];
		hi_hodo_ematch = new H1F[2];
		hi_hodo_tmatch = new H1F[2];
		hi_hodo_ematch_2D = new H2F[2];
		hi_hodo_tmatch_2D = new H2F[2];
		f_charge_landau = new F1D[2];

		for(int layer=0; layer < 2; layer++) {
			hi_hodo_eall[layer] = new H1F(String.format("hi_hodo_eall_l%d",layer+1), String.format("hi_hodo_eall_l%d",layer+1), 200, 0, 10);
			hi_hodo_eall[layer].setTitleX("E (MeV)");
			hi_hodo_eall[layer].setTitleY("Counts");
			hi_hodo_eall[layer].setFillColor(4);
			hi_hodo_ematch[layer] = new H1F(String.format("hi_hodo_ematch_l%d",layer+1), String.format("hi_hodo_ematch_l%d",layer+1),200, 0, 10);
			hi_hodo_ematch[layer].setTitleX(String.format("E (MeV)"));
			hi_hodo_ematch[layer].setTitleY(String.format("Counts"));
			hi_hodo_ematch[layer].setFillColor(3);
			f_charge_landau[layer] = new F1D(String.format("Landau_%d",layer+1),"[amp]*landau(x,[mean],[sigma])+[p0]+[p1]*x", 0.5*(layer+1), 10.0);
        		f_charge_landau[layer].setParameter(0,0.0);
                        f_charge_landau[layer].setParameter(1,0.0);
                        f_charge_landau[layer].setParameter(2,1.0);
                        f_charge_landau[layer].setParameter(3,0.0);
                        f_charge_landau[layer].setParameter(4,0.0);
                        f_charge_landau[layer].setOptStat(1111111);
                        f_charge_landau[layer].setLineWidth(2);
			hi_hodo_ematch_2D[layer] = new H2F(String.format("hi_hodo_ematch_2D_l%d",layer+1), String.format("hi_hodo_ematch_2D_l%d",layer+1), 100, 0, 10 , 118, 0, 118);
            		hi_hodo_ematch_2D[layer].setTitleX("E (MeV)");
            		hi_hodo_ematch_2D[layer].setTitleY("Tile");
			hi_hodo_tmatch[layer] = new H1F(String.format("hi_hodo_tmatch_l%d",layer+1),  String.format("hi_hodo_tmatch_l%d",layer+1), 100, -50, 50);
			hi_hodo_tmatch[layer].setTitleX(String.format("T-T_start (ns)"));
			hi_hodo_tmatch[layer].setTitleY(String.format("Counts"));
            		hi_hodo_tmatch[layer].setFillColor(3);
            		hi_hodo_tmatch_2D[layer] = new H2F(String.format("hi_hodo_tmatch_2D_l%d",layer+1),  String.format("hi_hodo_tmatch_2D_l%d",layer+1), 100, -50, 50, 118, 0, 118);
            		hi_hodo_tmatch_2D[layer].setTitleX("E (MeV)");
            		hi_hodo_tmatch_2D[layer].setTitleY("Tile");
		}

		//Calorimeter Histograms
		hi_cal_nclusters = new H1F("hi_cal_nclusters", "N. Clusters", "Counts", 5, 0, 5);
        	hi_cal_nclusters.setFillColor(44);
        	hi_cal_clsize = new H1F("hi_cal_clsize", "Cluster Size", "Counts", 25, 0, 25);
        	hi_cal_clsize.setFillColor(44);
        	hi_cal_clsize_ch = new H1F("hi_cal_clsize_ch", "Cluster Size", "Counts", 25, 0, 25);
        	hi_cal_clsize_ch.setFillColor(44);
        	hi_cal_clsize_en = new H2F("hi_cal_clsize_en", " ", 25, 0, 25, 100, 0, 12);
        	hi_cal_clsize_en.setTitleX("Cluster size");
        	hi_cal_clsize_en.setTitleY("E (GeV)");
        	hi_cal_e_all = new H1F("hi_cal_e_all", "E (GeV)", "Counts", 200, 0, 12);
        	hi_cal_e_all.setFillColor(4);
        	hi_cal_e_ch = new H1F("hi_cal_e_ch", "E (GeV)", "Counts", 200, 0, 12);
        	hi_cal_e_ch.setFillColor(2);
        	hi_cal_e_neu = new H1F("hi_cal_e_neu", "E (GeV)", "Counts", 200, 0, 12);
        	hi_cal_e_neu.setFillColor(3);
        	hi_cal_theta_ch = new H1F("hi_cal_theta_ch","#theta (deg)", "Counts", 100, 2,  6);
        	hi_cal_theta_ch.setFillColor(2);
        	hi_cal_phi_ch = new H1F("hi_cal_phi_ch", "#phi (deg)", "Counts", 100, -180,180);
        	hi_cal_phi_ch.setFillColor(2);
        	hi_cal_time_ch = new H1F("hi_cal_time_ch", "T-T_RF(ns)", "Counts", 100, -rfPeriod/2,rfPeriod/2);
        	hi_cal_time_ch.setFillColor(33);
        	hi_cal_time_cut_ch = new H1F("hi_cal_time_cut_ch", "T-T_RF(ns)", "Counts", 100, -rfPeriod/2,rfPeriod/2);
        	hi_cal_time_cut_ch.setFillColor(3);
        	ftime_ch = new F1D("ftime_ch", "[amp]*gaus(x,[mean],[sigma])", -1., 1.);
        	ftime_ch.setParameter(0, 0.0);
        	ftime_ch.setParameter(1, 0.0);
        	ftime_ch.setParameter(2, 2.0);
        	ftime_ch.setLineWidth(2);
        	ftime_ch.setOptStat("1111");
        	hi_cal_time_e_ch = new H2F("hi_cal_time_e_ch", "hi_cal_time_e_ch", 100, 0., 12., 100, -rfPeriod/2,rfPeriod/2);
        	hi_cal_time_e_ch.setTitleX("E (GeV)");
        	hi_cal_time_e_ch.setTitleY("T-T_RF (ns)");
        	hi_cal_time_theta_ch = new H2F("hi_cal_time_theta_ch", "hi_cal_time_theta_ch", 100, 2., 6., 100, -rfPeriod/2,rfPeriod/2);
        	hi_cal_time_theta_ch.setTitleX("#theta (deg)");
        	hi_cal_time_theta_ch.setTitleY("T-T_RF (ns)");
        	hi_cal_time_neu = new H1F("hi_cal_time_neu", "T-T_start(ns)", "Counts", 100, -2,2);
        	hi_cal_time_neu.setFillColor(44);
        	hi_cal_time_cut_neu = new H1F("hi_cal_time_cut_neu", "T-T_start(ns)", "Counts", 100, -2,2);
        	hi_cal_time_cut_neu.setFillColor(4);
        	ftime_neu = new F1D("ftime_neu", "[amp]*gaus(x,[mean],[sigma])", -1., 1.);
        	ftime_neu.setParameter(0, 0.0);
        	ftime_neu.setParameter(1, 0.0);
        	ftime_neu.setParameter(2, 2.0);
        	ftime_neu.setLineWidth(2);
        	ftime_neu.setOptStat("1111");
        	hi_cal_time_e_neu = new H2F("hi_cal_time_e_neu", "hi_cal_time_e_neu", 100, 0., 12., 100, -2,2);
        	hi_cal_time_e_neu.setTitleX("E (GeV)");
        	hi_cal_time_e_neu.setTitleY("T-T_start (ns)");
        	hi_cal_time_theta_neu = new H2F("hi_cal_time_theta_neu", "hi_cal_time_theta_neu", 100, 2., 6., 100, -2,2);
        	hi_cal_time_theta_neu.setTitleX("#theta (deg)");
        	hi_cal_time_theta_neu.setTitleY("T-T_start (ns)");

		//Pi0 Histograms
        	hpi0sum = new H1F("hpi0sum", 200,50., 250.);
        	hpi0sum.setTitleX("M (MeV)");
        	hpi0sum.setTitleY("Counts");
        	hpi0sum.setTitle("2#gamma invariant mass");
        	hpi0sum.setFillColor(3);
        	fpi0 = new F1D("fpi0", "[amp]*gaus(x,[mean],[sigma])+[p0]+[p1]*x", 80.,200.);
        	fpi0.setParameter(0, 0.0);
        	fpi0.setParameter(1, 140.0);
        	fpi0.setParameter(2, 2.0);
        	fpi0.setParameter(3, 0.0);
        	fpi0.setParameter(4, 0.0);
        	fpi0.setLineWidth(2);
        	fpi0.setOptStat("1111111");
        	hmassangle = new H2F("hmassangle", 100, 0., 300., 100, 0., 6.);
        	hmassangle.setTitleX("M (MeV)");
        	hmassangle.setTitleY("Angle (deg)");
        	hmassangle.setTitle("Angle vs. Mass");
		
	}

	public void fillFTHodo(DataBank HodoHits, DataBank HodoClusters) {
		for(int i=0; i<HodoHits.rows(); i++) {
			int hodoC = HodoHits.getShort("clusterID",i);
			int hodoS = HodoHits.getByte("sector",i);
			int hodoL = HodoHits.getByte("layer",i);
			int component = HodoHits.getShort("component",i);
			int tile = -1;
			switch (hodoS) {
				case 1:
                          	  tile = component + 0;
                          	  break;
				case 2:
                          	  tile = component + 9;
                          	  break;
				case 3:
                          	  tile = component + 29;
                          	  break;
				case 4:
                          	  tile = component + 38;
   				  break;
				case 5:
                          	  tile = component + 58;
                          	  break;
				case 6:
                          	  tile = component + 67;
                          	  break;
				case 7:
                          	  tile = component + 87;
                          	  break;
				case 8:
                          	  tile = component + 96;
                          	  break;
				default:
                          	  tile = -1;
                          	  break;
			}
			double hodoHitE = HodoHits.getFloat("energy",i);
			double hodoHitT = HodoHits.getFloat("time",i);
			double hodoHitX = HodoHits.getFloat("x",i);
			double hodoHitY = HodoHits.getFloat("y",i);
			double hodoHitZ = HodoHits.getFloat("z",i);
			double path = Math.sqrt(hodoHitX*hodoHitX+hodoHitY*hodoHitY+hodoHitZ*hodoHitZ);
			int   clusterId = HodoHits.getShort("clusterID",i);
			hi_hodo_eall[hodoL-1].fill(hodoHitE);

			for(int j=0; j<HodoClusters.rows(); j++) {
                        	if(clusterId==HodoClusters.getShort("id", j) && HodoClusters.getShort("size", j)>1) {
                            		hi_hodo_ematch[hodoL-1].fill(hodoHitE);
                            		hi_hodo_ematch_2D[hodoL-1].fill(hodoHitE,tile);
                            		if(startTime > -100) {
                                		hi_hodo_tmatch[hodoL-1].fill(hodoHitT-path/29.97-startTime);
                                		hi_hodo_tmatch_2D[hodoL-1].fill(hodoHitT-path/29.97-startTime,tile);
                            		}
                        	}
                    	}
                }
	}

	public void fillFTCalo(DataBank ftPart, DataBank CalClusters) {
		List<Particle> gammas = new ArrayList<Particle>();
		hi_cal_nclusters.fill(ftPart.rows());
		for (int loop = 0; loop < ftPart.rows(); loop++) {
                int    charge = ftPart.getByte("charge", loop);
                double energy = ftPart.getFloat("energy", loop);
                double   time = ftPart.getFloat("time", loop);
                double     cx = ftPart.getFloat("cx", loop);
                double     cy = ftPart.getFloat("cy", loop);
                double     cz = ftPart.getFloat("cz", loop);
                int     calID = ftPart.getShort("calID", loop);
                int    hodoID = ftPart.getShort("hodoID", loop);

                double energyR  = 0;
                int    size     = 0;
                double path     = 0;
                for(int i=0; i<CalClusters.rows(); i++) {
			if(calID == CalClusters.getShort("id", i)) {
                        energyR  = CalClusters.getFloat("recEnergy", i);
                        size     = CalClusters.getInt("size", i);
                        double x = CalClusters.getFloat("x", i);
                        double y = CalClusters.getFloat("y", i);
                        double z = CalClusters.getFloat("z", i);
                        path     = Math.sqrt(x*x+y*y+z*z);
                        time     = CalClusters.getFloat("time", i)-path/29.97;
                    }
                }
		hi_cal_clsize.fill(size);
                hi_cal_e_all.fill(energy);
                hi_cal_clsize_en.fill(size, energy);

		if (charge != 0) {
                    hi_cal_clsize_ch.fill(size);
                    hi_cal_e_ch.fill(energy);
                    hi_cal_theta_ch.fill(Math.toDegrees(Math.acos(cz)));
                    hi_cal_phi_ch.fill(Math.toDegrees(Math.atan2(cy,cx)));
                    if(rfTime!=-1000) {
			hi_cal_time_ch.fill((time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                        if(energy>2) hi_cal_time_cut_ch.fill((time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                        hi_cal_time_e_ch.fill(energy,(time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                        hi_cal_time_theta_ch.fill(Math.toDegrees(Math.acos(cz)),(time-rfTime+1000.5*rfPeriod)%rfPeriod-0.5*rfPeriod);
                    }
                }
		else {
                    Particle recParticle = new Particle(22, energy*cx, energy*cy, energy*cz, 0,0,0);
                    gammas.add(recParticle);
                    hi_cal_e_neu.fill(energy);
                    if(startTime!=-1000 && trigger==11) {
                        hi_cal_time_neu.fill(time-startTime);
                        if(energy>2) hi_cal_time_cut_neu.fill(time-startTime);
                        hi_cal_time_e_neu.fill(energy,time-startTime);
                        hi_cal_time_theta_neu.fill(Math.toDegrees(Math.acos(cz)),time-startTime);
                    }
                }
            }
		
	if(gammas.size()>=2) {
		for (int i1 = 0; i1 < gammas.size(); i1++) {
                	for (int i2 = i1 + 1; i2 < gammas.size(); i2++) {
                        	Particle partGamma1 = gammas.get(i1);
                        	Particle partGamma2 = gammas.get(i2);
                        	Particle partPi0 = new Particle();
                        	partPi0.copy(partGamma1);
                        	partPi0.combine(partGamma2, +1);
                        	double invmass = Math.sqrt(partPi0.mass2());
                        	double x = (partGamma1.p() - partGamma2.p()) / (partGamma1.p() + partGamma2.p());
                        	double angle = Math.toDegrees(Math.acos(partGamma1.cosTheta(partGamma2)));
                        	if(angle>1.5) hpi0sum.fill(invmass*1000);
                        	hmassangle.fill(invmass*1000, angle);
                    }
                }
            }

	}

        public void processEvent(DataEvent event) {
		
		DataBank recRun    = null;
        	DataBank recBankEB = null;
        	DataBank recEvenEB = null;
        	DataBank ftParticles = null;
        	DataBank ftCalClusters = null;
        	DataBank ftHodoClusters = null;
        	DataBank ftHodoHits = null;
        	if(event.hasBank("RUN::config"))            recRun      = event.getBank("RUN::config");
        	if(event.hasBank("REC::Particle"))          recBankEB   = event.getBank("REC::Particle");
        	if(event.hasBank("REC::Event"))             recEvenEB   = event.getBank("REC::Event");
        	if(event.hasBank("FT::particles"))          ftParticles = event.getBank("FT::particles");
        	if(event.hasBank("FTCAL::clusters"))      ftCalClusters = event.getBank("FTCAL::clusters");
        	if(event.hasBank("FTHODO::clusters"))    ftHodoClusters = event.getBank("FTHODO::clusters");
        	if(event.hasBank("FTHODO::hits"))            ftHodoHits = event.getBank("FTHODO::hits");

		// Decoding Trigger Bits
		if (event.hasBank("RUN::config")) {
			DataBank bank = event.getBank("RUN::config");
			long TriggerWord = bank.getLong("trigger", 0) & 0xFFFFFFFF;
			for (int i = 31; i >= 0; i--) {
				trigger_bits[i] = (TriggerWord & (1 << i)) != 0;
			}
		}

		//Get event start time
		if(recEvenEB!=null) {
			startTime = recEvenEB.getFloat("STTime", 0);
			rfTime    = recEvenEB.getFloat("RFTime", 0);
        	}

		//Get trigger particle
        	if(recBankEB!=null) trigger = recBankEB.getInt("pid", 0);

		//Main Processing
		if (ftParticles != null /*&& trigger_bits[25]*/) {
			if (ftHodoHits != null && ftHodoClusters != null) fillFTHodo(ftHodoHits, ftHodoClusters);
			if (ftCalClusters != null) fillFTCalo(ftParticles, ftCalClusters);
		} //End if ftParticle is not null
		

	}

	public void analyze() {
		//Fit hodoscope charge
        	for(int layer=0; layer < 2; layer++) {
            		initLandauFitPar(hi_hodo_ematch[layer], f_charge_landau[layer]);
            		DataFitter.fit(f_charge_landau[layer],hi_hodo_ematch[layer],"LRQ");
            		hi_hodo_ematch[layer].setFunction(null);
        	}
		//Fit calorimeter time
        	initTimeGaussFitPar(ftime_ch,hi_cal_time_cut_ch);
        	DataFitter.fit(ftime_ch,hi_cal_time_cut_ch,"LQ");
        	hi_cal_time_cut_ch.setFunction(null);
        	initTimeGaussFitPar(ftime_neu,hi_cal_time_cut_neu);
        	DataFitter.fit(ftime_neu,hi_cal_time_cut_neu,"LQ");
        	hi_cal_time_cut_neu.setFunction(null);
		//Fit pi0 mass
        	double hAmp  = hpi0sum.getBinContent(hpi0sum.getMaximumBin());
        	double hMean = hpi0sum.getAxis().getBinCenter(hpi0sum.getMaximumBin());
        	double hRMS  = 10; //ns
        	fpi0.setParameter(0, hAmp);
        	fpi0.setParLimits(0, hAmp*0.8, hAmp*1.2);
        	fpi0.setParameter(1, hMean);
        	fpi0.setParLimits(1, hMean-hRMS, hMean+hRMS);
        	DataFitter.fit(fpi0,hpi0sum,"LQ");
        	hpi0sum.setFunction(null);
	}

	private void initLandauFitPar(H1F hcharge, F1D fcharge) {
        	double hAmp  = hcharge.getBinContent(hcharge.getMaximumBin());
        	double hMean = hcharge.getAxis().getBinCenter(hcharge.getMaximumBin());
        	double hRMS  = hcharge.getRMS(); //ns
        	fcharge.setRange(fcharge.getRange().getMin(), hMean*2.0);
        	fcharge.setParameter(0, hAmp);
        	fcharge.setParLimits(0, 0.5*hAmp, 1.5*hAmp);
        	fcharge.setParameter(1, hMean);
        	fcharge.setParLimits(1, 0.8*hMean, 1.2*hMean);//Changed from 5-30        
        	fcharge.setParameter(2, 0.3);//Changed from 2
        	fcharge.setParLimits(2, 0.1, 1);//Changed from 0.5-10
        	fcharge.setParameter(3, 0.2*hAmp);
                fcharge.setParameter(4, -0.3);//Changed from -0.2
	}

	private void initTimeGaussFitPar(F1D ftime, H1F htime) {
        	double hAmp  = htime.getBinContent(htime.getMaximumBin());
        	double hMean = htime.getAxis().getBinCenter(htime.getMaximumBin());
        	double hRMS  = htime.getRMS(); //ns
        	double rangeMin = (hMean - (3*hRMS));
        	double rangeMax = (hMean + (3*hRMS));
        	double pm = hRMS*3;
        	ftime.setRange(rangeMin, rangeMax);
        	ftime.setParameter(0, hAmp);
        	ftime.setParLimits(0, hAmp*0.8, hAmp*1.2);
        	ftime.setParameter(1, hMean);
        	ftime.setParLimits(1, hMean-pm, hMean+(pm));
        	ftime.setParameter(2, 0.2);
        	ftime.setParLimits(2, 0.1*hRMS, 0.8*hRMS);
	}

        public void plot() {
                EmbeddedCanvas can_FT = new EmbeddedCanvas();
                can_FT.setSize(3000,5000);
                can_FT.divide(3,10);
                can_FT.setAxisTitleSize(18);
                can_FT.setAxisFontSize(18);
                can_FT.setTitleSize(18);
               	can_FT.cd(0);can_FT.draw(hi_hodo_eall[0]);can_FT.draw(hi_hodo_ematch[0],"same");can_FT.draw(f_charge_landau[0],"same");
		can_FT.cd(1);can_FT.getPad(1).getAxisZ().setLog(true);can_FT.draw(hi_hodo_ematch_2D[0]);
		can_FT.cd(3);can_FT.draw(hi_hodo_eall[1]);can_FT.draw(hi_hodo_ematch[1],"same");can_FT.draw(f_charge_landau[1],"same");
		can_FT.cd(4);can_FT.getPad(4).getAxisZ().setLog(true);can_FT.draw(hi_hodo_ematch_2D[1]);
		can_FT.cd(6);can_FT.draw(hi_hodo_tmatch[0]);
		can_FT.cd(7);can_FT.getPad(7).getAxisZ().setLog(true);can_FT.draw(hi_hodo_tmatch_2D[0]);
		can_FT.cd(9);can_FT.draw(hi_hodo_tmatch[1]);
		can_FT.cd(10);can_FT.getPad(10).getAxisZ().setLog(true);can_FT.draw(hi_hodo_tmatch_2D[1]);	
		can_FT.cd(12);can_FT.getPad(12).getAxisY().setLog(true);can_FT.draw(hi_cal_nclusters);
		can_FT.cd(13);can_FT.getPad(13).getAxisY().setLog(true);can_FT.draw(hi_cal_clsize);can_FT.draw(hi_cal_clsize_ch,"same");
		can_FT.cd(14);can_FT.getPad(14).getAxisZ().setLog(true);can_FT.draw(hi_cal_clsize_en);
		can_FT.cd(15);can_FT.draw(hi_cal_e_all);can_FT.draw(hi_cal_e_ch,"same");can_FT.draw(hi_cal_e_neu,"same");
		can_FT.cd(16);can_FT.draw(hi_cal_theta_ch);
		can_FT.cd(17);can_FT.draw(hi_cal_phi_ch);
		can_FT.cd(18);can_FT.draw(hi_cal_time_ch);can_FT.draw(hi_cal_time_cut_ch,"same");can_FT.draw(ftime_ch,"same");
		can_FT.cd(19);can_FT.draw(hi_cal_time_e_ch);
		can_FT.cd(20);can_FT.draw(hi_cal_time_theta_ch);
		can_FT.cd(21);can_FT.draw(hi_cal_time_neu);can_FT.draw(hi_cal_time_cut_neu,"same");can_FT.draw(ftime_neu,"same");
		can_FT.cd(22);can_FT.draw(hi_cal_time_e_neu);
		can_FT.cd(23);can_FT.draw(hi_cal_time_theta_neu);
		can_FT.cd(24);can_FT.draw(hpi0sum);can_FT.draw(fpi0,"same");
		can_FT.cd(25);can_FT.draw(hmassangle);
		

		if(runNum>0){
			if(!write_volatile)can_FT.save(String.format("plots"+runNum+"/FT.png"));
			if(write_volatile)can_FT.save(String.format("/volatile/clas12/rgb/spring19/plots"+runNum+"/FT.png"));
			System.out.println(String.format("saved plots"+runNum+"/FT.png"));
		}
		else{
			can_FT.save(String.format("plots/FT.png"));
			System.out.println(String.format("saved plots/FT.png"));
		}

	}
        public void write() {
                TDirectory dirout = new TDirectory();
                dirout.mkdir("/ft/");
                dirout.cd("/ft/");
                for(int s=0;s<2;s++){
                        dirout.addDataSet(hi_hodo_eall[s],hi_hodo_ematch[s],hi_hodo_ematch_2D[s],hi_hodo_tmatch[s],hi_hodo_tmatch_2D[s]);
		}
		dirout.addDataSet(hi_cal_nclusters,hi_cal_clsize,hi_cal_clsize_ch,hi_cal_clsize_en,hi_cal_e_ch,hi_cal_e_all,hi_cal_theta_ch,hi_cal_phi_ch,hi_cal_time_ch,hi_cal_time_cut_ch,hi_cal_time_e_ch);
		dirout.addDataSet(hi_cal_time_theta_ch,hi_cal_time_neu,hi_cal_time_cut_neu,hi_cal_time_e_neu,hi_cal_time_theta_neu,hpi0sum,hmassangle);
                if(write_volatile)if(runNum>0)dirout.writeFile("/volatile/clas12/rgb/spring19/plots"+runNum+"/out_FT_"+runNum+".hipo");
                
		if(!write_volatile){
			if(runNum>0)dirout.writeFile("plots"+runNum+"/out_FT_"+runNum+".hipo");
			else dirout.writeFile("plots/out_FT.hipo");
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
		FT ana = new FT(runNum,useTB,useVolatile);
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
			System.out.println(String.format(">>>>>>>>>>>>>>>> FT %s",runstrg));
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
				if(count%10000 == 0) System.out.println(count/1000 + "k events (this is FT on "+runstrg+") progress : " + progresscount+"/"+filetot);
			}
			reader.close();
		}
		System.out.println("Total : " + count + " events");
		ana.analyze();
		ana.plot();
		ana.write();
        }
}

