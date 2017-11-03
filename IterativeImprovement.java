package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.distribution.GeneTreeBrSpeciesNetDistribution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.likelihood.BeagleTreeLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.sitemodel.SiteModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.Frequencies;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.JukesCantor;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
//import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by doriswang on 3/12/17.
 */

//This class iteratively improve the likelihood of Gene tree. P(S|g)P(g|ST)

public class IterativeImprovement {
    private static String TREE_DIR = "/Users/doriswang/PhyloNet/Data/";
    private static String SEQ_DIR = "/Users/doriswang/PhyloNet/Data/";
    private static String RESULT_DIR = "/Users/doriswang/PhyloNet/Data/IIG/result/";
    private static String ST_DIR = "/Users/doriswang/PhyloNet/Data/";


    private Network<NetNodeInfo> INIT_ST; //cUnit
    private Network<NetNodeInfo> finalST;
    private Network<NetNodeInfo> bestEstimatedST;
    private String GLOBALBESTST;
    private String trueST;
    private String[] trueSTList;

    private List<String> trueGTS = new ArrayList<String>();
    private List<UltrametricTree> finalGTS = new ArrayList<UltrametricTree>();
    private List<String> GLOBALBESTGTS = new ArrayList<String>();

    private static List<Alignment> trueSeq = new ArrayList<Alignment>();
    private double GLOBALBESTLL;
    private IIGTSimulator simulator;
    private InferOperator operator;

    //TODO: PARA
    private static int ITERATION = 1000; //{1k, 5k ,10k}
    private static boolean ISREALDATA = false;
    private static boolean FULL_LL = true;
    private static String DATASET = "YEAST"; // "ROSE" "YEAST" "17"
    //private static double STREE_HEIGHT = 3.0; // "ROSE" -> 5.0 ELSE -> 3.0

    //"YEAST" : 7 taxa, 106 loci, 1 replication   || "ROSE" :  100 taxa, 25 loci, 100 replication, halfTheta 0.001
    //YEAST_halfTheta:  (DQ) 0.055   || (BEAST)   0.04      // (Luay) 7.568*10^-3 / 2 = 3.784*10^-3
    //17:2.5*10^(-9)
    private static double halfTheta = 0.04;//need to init in ALL classes ||  theta -> POP_SIZE

    private static Integer lociNum = 106;
    //#true gene trees
    private Integer trueGTSize = lociNum;
    //#internal gene tree size
    private static Integer refineGTSize = 50;
    private static double[] _scales = new double[]{1.0};
    private static int[] _seqLens = new int[]{2000};
    private static int taxaNum = 7;

    // LL(GTS|ST) for pseudo version || LL(Seq|GTS)*(GTS|ST) for classic version
    private double[] gtsLL = new double[ITERATION];// LL(Seq|GTS)*(GTS|ST)
    // LL(Seq|GTS)*(GTS|ST) for pseudo version || not used for classic version
    private double[] gtsHLL = new double[ITERATION];

    private double currentLL = 0.0;// LL(Seq|GTS)*(GTS|ST)

    private double currentHLL = 0.0;

    private double[] currentGTLL = new double[lociNum];
    private double[] currentHGTLL = new double[lociNum];

    private double[] gtDis = new double[ITERATION];
    private double[] stDis = new double[ITERATION];
    private double[][] stDis_Y = new double[ITERATION][4];
    private double maxLL;
    private int itNum;
    private long[] time = new long[2];
    private static String STRATEGY = "RANDOM"; //RANDOM IMPROVE NONE
    //private int data;


    IterativeImprovement() throws IOException, InterruptedException, ParseException{
        //TODO Simulation
        simulator = new IIGTSimulator(trueGTSize, _scales, _seqLens, halfTheta, ITERATION, "/Users/doriswang/PhyloNet/Data/IIG/");
        operator = new InferOperator("/Users/doriswang/PhyloNet/Data/IIG/","/Users/doriswang/PhyloNet/Data/IIG/output",0);
        String resultFolder = RESULT_DIR + ITERATION + "/" + taxaNum + "_" + lociNum + "/" + refineGTSize + "/";
        operator.isExitsPath(resultFolder);
        if (ISREALDATA) {
            if (DATASET.equals("YEAST")) {
                trueSeq = operator.loadYData("/Users/doriswang/PhyloNet/Data/IIG/seq/yeast_infer_network_from_multilocus_data_bayesian.nexus", 7, 106);
                //Nature; ((((Sbay,Spar)I5,Smik)I4,(Scer,Skud)I3)I2,(Scas,Sklu)I1);
                trueSTList = new String[4];
                //Nature & DQ(1)
                trueSTList[0] = "((((((Scer,Spar)I5,Smik)I4,Skud)I3,Sbay)I2,Scas)I1,Sklu);";
                //*beast 1
                trueSTList[1] = "((((Sbay,Spar)I5,Smik)I4,(Scer,Skud)I3)I2,(Scas,Sklu)I1)";
                //*beast 2
                trueSTList[2] = "(((((Sbay,Spar)I5,Smik)I4,(Skud,Scer)I3)I2,Scas)I1,Sklu);";
                //DQ(2)
                trueSTList[3] = "(((((Scer,Spar)I5,Smik)I4,(Skud,Sbay)I3)I2,Scas)I1,Sklu);";
                //stDis_Y = new double[6][];
            }
            if (DATASET.equals(("ROSE"))) {
                //TODO 6.5 random repNumber(0-9)
                //int repNum = random
                operator.loadRData("/Users/doriswang/PhyloNet/Data/ROSE/seq/", 0, 25, trueSeq, trueGTS);
                trueST = "((((SIE:32.768400,(SHH:55.161100,SD:92.887200):0.017200):147.054800,(SGF:313.737600,((SDD:136.012900,SGJ:185.723500):142.756300,((SDA:193.043900,(SIC:13.622900,SGB:17.998700):217.978200):41.139300,((SEF:123.695200,SEE:74.433100):8.837600,((SGE:40.694700,SCA:32.719600):3.969200,SHF:55.601100):51.256000):72.577700):104.657600):0.167100):24.663100):436.564100,((SFB:301.990200,(SC:241.140400,(SBB:80.283200,SBA:154.224300):48.452400):270.521400):89.901500,(SII:252.893300,((SBAA:529.313000,((SJA:94.094300,(SIA:87.710400,(SJC:58.506500,SJJ:71.226900):23.032900):146.525600):120.254700,((SEI:65.299400,SJD:103.739600):588.950600,(((((SCJ:549.947900,((((SDE:66.548000,SCG:37.359000):40.193500,(SDF:154.793700,SDC:145.448000):41.021400):1.755400,(((SHC:4.300900,SBE:6.576000):10.713400,SFJ:5.852000):2.525500,SIG:19.862200):78.497800):299.186900,(SDH:171.889200,((SFC:51.965800,SCD:100.723900):16.448000,((SCB:105.204000,SHI:75.456400):15.978200,SB:51.330300):32.436700):10.222800):236.314800):23.223600):12.447700,((SI:251.074600,(SJI:48.661600,(SEH:18.160700,SEA:13.852200):70.756800):194.742300):109.815700,SFD:291.136400):6.360800):14.906900,(SIB:13.180200,SGC:30.112900):270.164800):23.633800,(((((SBH:101.064700,((SEJ:22.118200,SFA:35.905100):21.071700,(SJB:16.228500,SCF:39.466000):10.098200):16.527900):30.307100,(((SEC:5.166000,SCH:7.083200):1.255200,SG:5.779700):6.457000,SEB:14.848900):44.347700):1.345100,SE:50.936000):78.367000,SJG:100.327200):180.359300,((((SCE:8.226400,SJE:9.550700):18.353100,SGH:28.108000):33.063300,(SIH:58.494700,SBG:42.460000):10.473200):139.302200,(SHE:199.469000,(SJF:43.082500,(SID:21.613800,SBF:18.909800):18.984000):89.915000):239.073300):54.232800):8.426900):8.148800,SFE:645.334900):1.739600):24.065000):6.917300):1.635300,((SDJ:19.966600,SGI:17.856800):286.158100,((SFG:138.272400,((SJ:9.085400,SBC:8.840900):15.029600,SBJ:13.470600):135.651700):75.775700,(((SBD:79.405100,SFH:38.904400):267.784100,SIJ:138.749700):26.762800,((SIF:90.883100,SHB:75.700800):373.834800,(SCI:75.460600,(SDG:99.975900,SDB:41.332600):64.678700):216.202800):9.818800):31.403700):29.469100):12.027700):20.852000):31.251300):156.810400):15.134200,(((SHA:184.003800,((((SHJ:1.748200,SED:3.337200):15.234200,SFF:22.713900):76.115700,SFI:135.936200):217.175600,SEG:261.659500):29.998500):153.458800,(SCC:187.929900,(SHG:4.238200,SGD:7.222000):93.105300):576.901900):2.722200,((SJH:24.757100,SGG:7.702900):238.807200,(SF:267.369500,((SH:76.450300,SBI:122.608800):156.785700,(SDI:26.352100,(SGA:22.043200,SHD:10.069600):23.156200):363.368300):120.732000):21.289200):139.951700):321.529500);";
            }

            if (DATASET.equals(("17"))) {
                //TODO 6.5 random repNumber(0-9)
                //int repNum = random

                trueST = operator.load17Data("/Users/doriswang/PhyloNet/Data/17-taxon/", 1000, 32, trueSeq, trueGTS);
            }
        }
     else {
            List<String> trees = simulator.simulateData();

            trueST = trees.get(0);
            trueGTS = trees.subList(1, trees.size());
//            for (int i = 0; i < lociNum; i++) {
//                Alignment aln = operator.loadLocus(i, _seqLens[0], taxaNum);
//                trueSeq.add(aln);
//            }
        }
    }
//    public List<Alignment> loadData(int stNum) throws IOException{
//        String gtFileName =  "/Users/doriswang/PhyloNet/Data/17-taxon/32loci/" + "Rep" + stNum + "gtrees";
//        String seqFileName = "/Users/doriswang/PhyloNet/Data/17-taxon/32loci/seq/" + "Rep" + stNum + "gseqs";
//        String stFileName =  "/Users/doriswang/PhyloNet/Data/17-taxon/st/" + "Rep" + stNum + "stree";
//
//        BufferedReader stReader = new BufferedReader(new FileReader(stFileName));
//        String stString = stReader.readLine();
//        Network truest = Networks.readNetwork(stString);
//        this.trueST = stString;
//
//        BufferedReader gtReader = new BufferedReader(new FileReader(gtFileName));
//        List<Tree> trees = new ArrayList<Tree>();
//        for(int i = 0; i<lociNum; i++){
//            Tree tempGT = Trees.readTree(gtReader.readLine());
//            trees.add(tempGT);
//        }
//        this.trueGTS =  operator.getUTreeByTree(trees);
//
//        return operator.loadAlnSQ(seqFileName,_seqLens[0], taxaNum,lociNum);
//
//    }

    public Network<NetNodeInfo> initST(List<Alignment> aln) throws IOException, ParseException {
        List tempGTs = operator.simGTSByUPGMA(aln); // Time: tao
        if (DATASET.equals("YEAST")) {
            for (int i = 0; i < tempGTs.size(); i++) {
                trueGTS.add(tempGTs.get(i).toString());
            }
        }
        operator.scaleGTS(tempGTs, 1 / halfTheta, false); //Time: c Unit
//        double height = operator.getNetHeight(tempST1);
        //Network<NetNodeInfo> tempST = operator.inferSTByGLASS(tempGTs, halfTheta);
        Network<NetNodeInfo> tempST = operator.inferSTByAST(tempGTs, halfTheta);

        operator.scaleSTNet(tempST, halfTheta, true);
        System.out.println("temp_ST-" + itNum + "-:" + tempST.toString());
        //double tHeight = operator.getNetHeight(Networks.readNetwork(trueST));
        currentGTLL = new double[lociNum];
        //double ll = computeFullLL(tempST, tempGTs);
        computeFullLL(tempST, tempGTs);
        gtsLL[0] = this.currentLL;
        gtsHLL[0] = this.currentHLL;
                //currentLL = gtsLL[0];
        maxLL = gtsLL[0];
        if(!FULL_LL)
            maxLL = gtsHLL[0];
        INIT_ST = tempST;
        System.out.println("Init_ST-" + itNum + "-:" + INIT_ST.toString());
        finalGTS = tempGTs;
        bestEstimatedST = INIT_ST;
        GLOBALBESTST = INIT_ST.toString();
        Iterator it = tempGTs.iterator();
        GLOBALBESTGTS = new ArrayList<String>();
        while (it.hasNext()) {
            GLOBALBESTGTS.add(it.next().toString());
        }
        GLOBALBESTLL = maxLL;
        Tree temp = Trees.readTree(GLOBALBESTST);
        //temp.getNode("A");
        if(DATASET.equals("YEAST")) {
            //stDis_Y = new double[ITERATION][6];
            for (int i = 0; i < 4; i++) {
                stDis_Y[0][i] = operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueSTList[i]));

            }
        }
        else
            stDis[0] = operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));

        double thisD = 0.0;
        for (int k = 0; k < lociNum; k++) {
            thisD += operator.getDistance((Tree) finalGTS.get(k).getTree(), Trees.readTree(trueGTS.get(k)));
        }
        gtDis[0] = thisD / lociNum;
        return INIT_ST;
    }


    //Input: alignments & t    Output: GTs and STs
    public List<UltrametricTree> iigt(List<Alignment> aln, int t) throws IOException, ParseException, InterruptedException {
        //CHANGE: PARA
        itNum = 1;
        ITERATION = t;
        long start = System.currentTimeMillis();

        //Initialization
        Network<NetNodeInfo> tempST = initST(aln);
        for (int i = 1; i < t; i++) {
            System.out.println("\n" + "#" + itNum + " iteration start! ");
            refineGeneSet(tempST);
            itNum++;
            //REMAINING: scale by height
//            if (itNum % 50 == 0) {
//                System.out.println("Before--50--Rescale : " + tempST.toString());
//                double currentHeight = 0.0;
//                NetNode leaf = (NetNode) tempST.getLeaves().iterator().next();
//                while (leaf != tempST.getRoot()) {
//                    currentHeight += leaf.getParentDistance((NetNode) leaf.getParents().iterator().next());
//                    leaf = (NetNode) leaf.getParents().iterator().next();
//                }
//                double ratio = STREE_HEIGHT / currentHeight;
//                operator.scaleSTNet(tempST, ratio, false);
//                System.out.println("After--50--Rescale : " + tempST.toString());
//            }
        }
        long end = System.currentTimeMillis();
        long costtime = end - start;
        String resultFolder = RESULT_DIR + ITERATION + "/" + taxaNum + "_" + lociNum + "/" + refineGTSize + "/";
        operator.isExitsPath(resultFolder);
        BufferedWriter llOut1 = new BufferedWriter(new FileWriter(resultFolder + "RunningTime.txt"));
        llOut1.write("Running time:" + String.valueOf(costtime) + "\n");
        llOut1.write("Running time for p(GTS|ST):" + String.valueOf(time[0]) + "\n");
        llOut1.write("Running time for p(Seq|GTS)*p(GTS|ST):" + String.valueOf(time[1]) + "\n");

        System.out.println("------Running time: " + costtime);
        //operator.isExitsPath(resultFolder);
        llOut1.flush();
        llOut1.close();
        BufferedWriter llOut = new BufferedWriter(new FileWriter(resultFolder + "Likelihood.txt"));
        for (int k = 0; k < ITERATION; k++) {
            llOut.write(k + " " + String.valueOf(gtsHLL[k]) + "\n");
        }
        llOut.flush();
        llOut.close();
        BufferedWriter llPOut = new BufferedWriter(new FileWriter(resultFolder + "FullLikelihood.txt"));
        for (int k = 0; k < ITERATION; k++) {
            llPOut.write(k + " " + String.valueOf(gtsLL[k]) + "\n");
        }
        llPOut.flush();
        llPOut.close();

        BufferedWriter bw1 = new BufferedWriter(new FileWriter(resultFolder + "BestST" + ".txt"));
        bw1.write(bestEstimatedST.toString() + "\n");
        bw1.flush();
        bw1.close();

        for (int tid = 0; tid < finalGTS.size(); tid++) {
            String gt = finalGTS.get(tid).toString();
            try {
                String fileName = resultFolder + "BestGT_locus" + tid + ".tree";
                BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                out.write(gt + "\n");
                out.flush();
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        //String g_resultFolder = RESULT_DIR + ITERATION + "/global/";
        BufferedWriter bw2 = new BufferedWriter(new FileWriter(resultFolder + "GlobalBestST" + ".txt"));
        bw2.write(GLOBALBESTST + "\n");
        bw2.write("Best Likelihood value : " + String.valueOf(GLOBALBESTLL) + "\n");
        bw2.write("Number of loci : " + lociNum + "\n");
        bw2.write("Number of refine gene tree set :  " + refineGTSize + "\n");
        bw2.write("Number of iteration :  " + ITERATION + "\n");
        bw2.write("SeqLens :  " + _seqLens[0] + "\n");
        bw2.write("Is real data :  " + ISREALDATA + "\n");
        bw2.write("Is full likelihood  " + FULL_LL + "\n");
        bw2.write("HalfTheta  " + halfTheta + "\n");
        bw2.flush();
        bw2.close();

        for (int tid = 0; tid < GLOBALBESTGTS.size(); tid++) {
            String gt = GLOBALBESTGTS.get(tid);
            try {
                String fileName = resultFolder + "GlobalBestGT_locus" + tid + ".tree";
                BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                out.write(gt + "\n");
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        BufferedWriter ldOut = new BufferedWriter(new FileWriter(resultFolder + "distance.txt"));
        //6.10 TODO multi-ST distance
        double totalD = 0.0;


        if(DATASET.equals("YEAST")){
            for(int i=0;i<4;i++){
                    ldOut.write("RF of ST " + String.valueOf(i) + " :" + String.valueOf(operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueSTList[i]))) + "\n");

//                ldOut.flush();
//                ldOut.close();
            }
        }
        else {
            String stDistance = "RF of ST : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
            ldOut.write(stDistance + "\n");
        }

            String gtDistance = "RF of GT_";
            double thisD = 0.0;
        for (int k = 0; k < lociNum; k++) {
            thisD = operator.getDistance(Trees.readTree(GLOBALBESTGTS.get(k)), Trees.readTree(trueGTS.get(k)));
            ldOut.write(gtDistance + k + " :" + thisD + "\n");
            totalD += thisD;
        }
        ldOut.write("Total distance of GTS is " + totalD + "\n");
        ldOut.write("Average distance of GTS is " + totalD / trueGTS.size() + "\n");
        BufferedWriter gtD = new BufferedWriter(new FileWriter(resultFolder + "GTdistance.txt"));
        for (int gtn = 0; gtn < ITERATION; gtn++) {
            gtD.write(String.valueOf(gtn) + ":" + String.valueOf(gtDis[gtn]) + "\n");
        }
        gtD.flush();
        gtD.close();
        ldOut.flush();
        ldOut.close();
        //6.10 TODO multi-ST distance
        if(DATASET.equals("YEAST")){
            for(int i=0;i<4;i++){
                BufferedWriter stD = new BufferedWriter(new FileWriter(resultFolder + "STdistance_" + String.valueOf(i) + ".txt"));
                for (int stn = 0; stn < ITERATION; stn++) {
                    stD.write(String.valueOf(stn) + ":" + String.valueOf(stDis_Y[stn][i]) + "\n");
                }
                stD.flush();
                stD.close();
            }
        }
        else {
            BufferedWriter stD = new BufferedWriter(new FileWriter(resultFolder + "STdistance.txt"));
            for (int stn = 0; stn < ITERATION; stn++) {
                stD.write(String.valueOf(stn) + ":" + String.valueOf(stDis[stn]) + "\n");
            }
            stD.flush();
            stD.close();
        }
        return finalGTS;
    }

    public List<UltrametricTree> refineGeneSet(Network<NetNodeInfo> tempST) throws ParseException, IOException {

        List<UltrametricTree> currentGTS = simGTSByMS(tempST, halfTheta, refineGTSize); //operator.MS
        List<UltrametricTree> gts = new ArrayList<>();
        currentGTLL = new double[lociNum];
        double ll = 0.0;
        double hLL = 0.0;
        List<Double> gProST = new ArrayList<>();
        if(itNum%5==0){
            //java.lang.OutOfMemoryError: Java heap space for 100 node case
            gProST = getGTSLLByST(currentGTS, tempST);
        }
        else
            gProST = getGTSLLByST(currentGTS, tempST);
        for (int i = 0; i < lociNum; i++) {
                UltrametricTree bestGTi = getBestGT(i, currentGTS, gProST);
                gts.add(bestGTi);
                ll += currentGTLL[i];
                hLL += currentHGTLL[i];
            }
        currentLL = ll;
        currentHLL = hLL;

        tempST = operator.inferSTByGLASS(gts, halfTheta);
        //tempST = operator.inferSTByAST(gts, halfTheta);
        boolean move = false;
        if(!FULL_LL)
            ll=hLL;
        if (STRATEGY == "NONE") {
            move = true;
        } else if (STRATEGY == "IMPROVE") {
            if (ll > maxLL) {
                move = true;
                GLOBALBESTST = tempST.toString();
                Iterator it = gts.iterator();
                GLOBALBESTGTS = new ArrayList<String>();
                while (it.hasNext()) {
                    GLOBALBESTGTS.add(it.next().toString());
                    GLOBALBESTLL = ll;
                }
            } else
                move = false;
        } else if (STRATEGY == "RANDOM") {
            if (ll > maxLL) {
                move = true;
                if (ll > GLOBALBESTLL) {
                    GLOBALBESTST = tempST.toString();
                    Iterator it = gts.iterator();
                    GLOBALBESTGTS = new ArrayList<String>();
                    while (it.hasNext()) {
                        GLOBALBESTGTS.add(it.next().toString());
                    }
                    GLOBALBESTLL = ll;
                }
            } else {
                double random = Math.random();
                if (random < 0.975)
                    move = false;
                else {
                    System.out.println("--------Random-Accept: #Iteration " + itNum + " : " + ll);
                    move = true;
                }
            }
        }
        if (move == true) {
            //6.11 TODO: scale by ll or half_LL
            if(FULL_LL)
                ll = scaleAll(gts, tempST, ll);
            else{
                ll = scaleAllHalf(gts, tempST, ll);
            }
//            gtsLL[itNum] = ll;
//            maxLL = ll;
//            finalGTS = gts;
//            bestEstimatedST = tempST;
            if(DATASET.equals("YEAST")) {
                for (int i = 0; i < 4; i++) {
                    stDis_Y[itNum][i] = operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueSTList[i]));
                    System.out.println("ST Distance of # " + itNum + "-" + i + " is : " + stDis_Y[itNum][i]);
                }

            }
            else {
                stDis[itNum] = operator.getDistance(Trees.readTree(tempST.toString()), Trees.readTree(trueST));
                System.out.println("ST Distance of # " + itNum + " is : " + stDis[itNum]);
            }
            double thisD = 0.0;
            for (int k = 0; k < lociNum; k++) {
                if(k==106)
                    System.out.print("l");
                thisD += operator.getDistance(gts.get(k).getTree(), Trees.readTree(trueGTS.get(k)));
            }
            gtDis[itNum] = thisD / lociNum;
            System.out.println("GTS Distance of # " + itNum + " is : " + gtDis[itNum]);
            System.out.println("**********Accept ST  : " + tempST.toString());

        } else {
            System.out.println("----------Reject ST  : " + tempST.toString());
            gtsLL[itNum] = gtsLL[itNum-1];
            gtsHLL[itNum] = gtsHLL[itNum-1];
            //6.11 TODO: scale by ll or half_LL
//            gtsLL[itNum] = maxLL;
//            if (!ISFULL_LL)
//                if(FULL_LL_NEED)
//                    gtsFLL[itNum] = fLL;
            if(DATASET.equals("YEAST")) {
                for (int i = 0; i < 4; i++) {
                    stDis_Y[itNum][i] = stDis_Y[itNum-1][i];

                }
            }
            else {
                stDis[itNum] = stDis[itNum-1];
                System.out.println("ST Distance of # " + itNum + " is : " + stDis[itNum]);
            }

            gtDis[itNum] = gtDis[itNum-1];
            System.out.println("GTS Distance of # " + itNum + " is : " + gtDis[itNum]);

        }

        System.out.println(gtsLL[itNum] + " :   : " + gtsHLL[itNum]);
        return finalGTS;
    }



//    //CORE IDEA: return one best gt for alignment_i
//    public UltrametricTree getBestGTYun(int locus, Network<NetNodeInfo> currentST, List<UltrametricTree> currentGeneSet) {
//        double bestLL = 0.0;
//        Alignment aln = trueSeq.get(locus);
//        UltrametricTree bestGeneTree = null;
//        for(int i = 0; i < currentGeneSet.size(); i++){
//            UltrametricTree gt = currentGeneSet.get(i);
//            double p1 = getGTLLByST(gt, currentST);
//            double p2 = getSeqLLByGT(gt, aln);
//            double p = p1*p2;
//            if(p >= bestLL) {
//                bestGeneTree = gt;
//                bestLL = p;
//            }
//        }
//        currentGTLL[locus] = bestLL;
//        return bestGeneTree;
//    }
//
//    //CORE IDEA: return one best gt for alignment_i
//    public UltrametricTree getBestGT(int locus, Network<NetNodeInfo> currentST, List<UltrametricTree> currentGeneSet) {
//        double bestLL = 0.0;
//        Alignment aln = trueSeq.get(locus);
//        UltrametricTree bestGeneTree = null;
//        for(int i = 0; i < currentGeneSet.size(); i++){
//            UltrametricTree gt = currentGeneSet.get(i);
//            double p1 = getGTLLByST(gt, currentST);
//            double p2 = getSeqLLByGT(gt, aln);
//            double p = p1*p2;
//            if(p >= bestLL) {
//                bestGeneTree = gt;
//                bestLL = p;
//            }
//        }
//        currentGTLL[locus] = bestLL;
//        return bestGeneTree;
//    }

    //CORE IDEA: return one best gt for alignment_i
    public UltrametricTree getBestGT(int locus, List<UltrametricTree> currentGeneSet, List<Double> gProST) throws IOException{
        double bestLL = Double.NEGATIVE_INFINITY;
        double bestHLL = Double.NEGATIVE_INFINITY;

        Alignment aln = trueSeq.get(locus);
        UltrametricTree bestGeneTree = currentGeneSet.get(0);
        for(int i = 0; i < currentGeneSet.size(); i++){
            List<UltrametricTree> uList1 = new ArrayList<UltrametricTree>();
            uList1.add(currentGeneSet.get(i));
            UltrametricTree gt = operator.scaleGTS(uList1,halfTheta,false).iterator().next();

            double p2 = getSeqLLByGT(gt, aln);
            if(p2>0) {
                p2 = getSeqLLByGT(gt, aln);
            }

            if(p2>0) {
                System.out.println("p(Seq|GT):" + p2 + "  || Locus: " + locus);
                continue;
            }
            double p = gProST.get(i) + p2;
            if(FULL_LL) {
                if (p >= bestLL) {
                    bestGeneTree = gt;
                    bestLL = p;
                    currentGTLL[locus] = bestLL;
                    currentHGTLL[locus] = p2;

                }
            }
            else{
                if (p2 >= bestLL) {
                    bestGeneTree = gt;
                    bestLL = p2;
                    currentHGTLL[locus] = bestLL;
                    currentGTLL[locus] = p;
                }
            }
        }
        //currentGTLL[locus] = bestLL;

        operator.scaleGTS(currentGeneSet,halfTheta,true).iterator().next();
//        if(badP2) {
//            System.out.println("Bad p2: best GT is " + bestGeneTree.toString());
//            System.out.println("Bad p2: init GT is " + initGTS.get(locus).toString());
//            System.out.println("-----Bad P End ------------------------------");
//        }
        return bestGeneTree;
    }

//    //CORE IDEA: return one best gt for alignment_i
//    public UltrametricTree getBestGTBYSeq(int locus, List<UltrametricTree> currentGeneSet) throws IOException{
//        double bestLL = Double.NEGATIVE_INFINITY;
//        boolean badP2 = false;
//        Alignment aln = trueSeq.get(locus);
//        UltrametricTree bestGeneTree = currentGeneSet.get(0);
//        for(int i = 0; i < currentGeneSet.size(); i++){
//            List<UltrametricTree> uList1 = new ArrayList<UltrametricTree>();
//            uList1.add(currentGeneSet.get(i));
//            UltrametricTree gt = operator.scaleGTS(uList1,halfTheta,false).iterator().next();
//            //TODO 6.3: delete p2
//            double p2 = getSeqLLByGT(gt, aln);
//            if(p2>0) {
//                p2 = getSeqLLByGT(gt, aln);
//            }
//            double p = (double)gProST.get(i) + p2;
//            if(p>0) {
//                //System.out.println("p:" + p + "  || Locus: " + locus);
//                continue;
//            }
//            if(p >= bestLL) {
//                if(p2>0)
//                    continue;
//                bestGeneTree = gt;
//                bestLL = p;
//            }
//        }
//        currentGTLL[locus] = bestLL;
//        operator.scaleGTS(currentGeneSet,halfTheta,true).iterator().next();
////        if(badP2) {
////            System.out.println("Bad p2: best GT is " + bestGeneTree.toString());
////            System.out.println("Bad p2: init GT is " + initGTS.get(locus).toString());
////            System.out.println("-----Bad P End ------------------------------");
////        }
//        return bestGeneTree;
//    }

    //CORE IDEA: return one best gt for alignment_i
//    public UltrametricTree getBestGTBySTYF(int locus, List<Double> gProST, List<UltrametricTree> currentGeneSet) throws IOException{
//        double bestLL = Double.NEGATIVE_INFINITY;
//        //boolean badP2 = false;
//        Alignment aln = trueSeq.get(locus);
//        int bestGTNum = 0;
//        UltrametricTree bestGeneTree = currentGeneSet.get(0);
//        for(int i = 0; i < currentGeneSet.size(); i++){
//            List<UltrametricTree> uList1 = new ArrayList<UltrametricTree>();
//            uList1.add(currentGeneSet.get(i));
//            UltrametricTree gt = operator.scaleGTS(uList1,halfTheta,false).iterator().next();
//            //TODO 6.3: delete p2
//            double p2 = getSeqLLByGT(gt, aln);
//            if(p2>0) {
//                p2 = getSeqLLByGT(gt, aln);
//            }
//            if(p2>0) {
//                continue;
//            }
//            double p = (double)gProST.get(i) + p2;
//            if(p>0) {
//                //System.out.println("p:" + p + "  || Locus: " + locus);
//                continue;
//            }
//            if(p >= bestLL) {
//                bestGeneTree = gt;
//                bestLL = p;
//            }
//        }
////        for(int i = 0; i < currentGeneSet.size(); i++){
////            List<UltrametricTree> uList1 = new ArrayList<UltrametricTree>();
////            uList1.add(currentGeneSet.get(i));
////            UltrametricTree gt = operator.scaleGTS(uList1,halfTheta,false).iterator().next();
////            //TODO 6.3: delete p2
////            double p = (double)gProST.get(i) + getSeqLLByGT(gt, aln);
////            if(p>0) {
////                //System.out.println("p:" + p + "  || Locus: " + locus);
////                continue;
////            }
////            if(p >= bestLL) {
////                bestLL = p;
////                bestGTNum = i;
////            }
////        }
//        currentFGTLL[locus] = bestLL;
//        operator.scaleGTS(currentGeneSet,halfTheta,true).iterator().next();
//
//        currentGTLL[locus] = gProST.get(bestGTNum);
//        return bestGeneTree;
//    }

    //If Full_LL return fullLL   Else Return half
    public double computeFullLL(Network<NetNodeInfo> currentST, List<UltrametricTree> currentGTS) {
        double ll = 0.0;
        double hll = 0.0;
        List<Double> gProST = getGTSLLBySTYF(currentGTS, currentST);
        for (int locus = 0; locus < trueSeq.size(); locus++) {
            Alignment aln = trueSeq.get(locus);
            UltrametricTree tempGT = currentGTS.get(locus);
            Iterator itNode1 = tempGT.getTree().getNodes().iterator();
            while (itNode1.hasNext()) {
                STINode n = (STINode) itNode1.next();
                //n.getParentDistance()
                if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(n.getParentDistance() * halfTheta); // Coalescent_Unit = Tao /  theta/2
                }
            }
            double p2 = getSeqLLByGT(tempGT, aln);
            double p = gProST.get(locus) + p2;

            if (p2 > 0) {
                System.out.println("p2: " + p2);
                continue;
            }
            if (p > 0) {
                System.out.println("p:" + p);
                continue;
            }
            ll += p;
            hll += p2;
            Iterator itNode2 = tempGT.getNodes().iterator();
            while (itNode2.hasNext()) {
                STINode n = (STINode) itNode2.next();
                //n.getParentDistance()
                if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(n.getParentDistance() / halfTheta); // Coalescent_Unit = Tao /  theta/2
                }
            }
        }
        this.currentLL = ll;
        this.currentHLL = hll;
        return ll;
    }

    //If Full_LL return fullLL   Else Return half
    public double computeHalfLL(Network<NetNodeInfo> currentST, List<UltrametricTree> currentGTS) {
        double ll = 0.0;
        double hll = 0.0;
        List<Double> gProST = getGTSLLBySTYF(currentGTS, currentST);
        for (int locus = 0; locus < trueSeq.size(); locus++) {
            Alignment aln = trueSeq.get(locus);
            UltrametricTree tempGT = currentGTS.get(locus);
            Iterator itNode1 = tempGT.getTree().getNodes().iterator();
            while (itNode1.hasNext()) {
                STINode n = (STINode) itNode1.next();
                //n.getParentDistance()
                if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(n.getParentDistance() * halfTheta); // Coalescent_Unit = Tao /  theta/2
                }
            }
            double p2 = getSeqLLByGT(tempGT, aln);
            double p = gProST.get(locus) + p2;

            if (p2 > 0) {
                System.out.println("p2: " + p2);
                continue;
            }
            ll += p;
            hll += p2;
            Iterator itNode2 = tempGT.getNodes().iterator();
            while (itNode2.hasNext()) {
                STINode n = (STINode) itNode2.next();
                //n.getParentDistance()
                if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(n.getParentDistance() / halfTheta); // Coalescent_Unit = Tao /  theta/2
                }
            }
        }
        this.currentLL = ll;
        this.currentHLL = hll;
        return hll;
    }


    public double computeGTLLByST(Network<NetNodeInfo> currentST, List<UltrametricTree> currentGTS){
        double ll = 0.0;
        //6.3 TODO compute full LL
        List<Double> gProST = getGTSLLBySTYF(currentGTS, currentST);
        for(int locus = 0 ; locus < trueSeq.size() ; locus++) {
            ll += gProST.get(locus);
        }
        return ll;
    }
//
//    public Network<NetNodeInfo> scaleST(Network<NetNodeInfo> currentST){
//        //for(each node)
//        double currentHeight = 0.0;
//        NetNode leaf = (NetNode)currentST.getLeaves().iterator().next();
//        while (leaf!=currentST.getRoot()){
//            currentHeight+=leaf.getParentDistance((NetNode) leaf.getParents().iterator().next());
//            leaf = (NetNode) leaf.getParents().iterator().next();
//        }
//        double ratio = STREE_HEIGHT/currentHeight;
//        Iterator itNode = currentST.getTreeNodes().iterator();
//        while (itNode.hasNext()) {
//            NetNode n = (NetNode) itNode.next();
//            if(n == (NetNode) currentST.getRoot())
//                continue;
//            //n.getParentDistance()
//            NetNode p = (NetNode) n.getParents().iterator().next();
//            if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
//                n.setParentDistance(p,n.getParentDistance(p) * ratio); // Coalescent_Unit = Tao /  theta/2
//            }
//        }
//
//        return currentST;
//    }

//
    public Network<NetNodeInfo> scaleSTNode(Network<NetNodeInfo> currentST){
        double currentHeight = 0.0;
        NetNode leaf = (NetNode)currentST.getLeaves().iterator().next();
        while (leaf!=currentST.getRoot()){
            currentHeight+=leaf.getParentDistance((NetNode) leaf.getParents().iterator().next());
            leaf = (NetNode) leaf.getParents().iterator().next();
        }
        Iterator itNode = currentST.getTreeNodes().iterator();
        while (itNode.hasNext()) {
            NetNode n = (NetNode) itNode.next();
            if(n == (NetNode) currentST.getRoot())
                continue;
            double heightN = 0.0;
            while(!n.isLeaf()) {
                if (n.getChildren().iterator().hasNext()) {
                    NetNode child = (NetNode) n.getChildren().iterator().next();
                    heightN += child.getParentDistance(n);
                    n = child;
                }
            }

            NetNode p = (NetNode) n.getParents().iterator().next();
            double ratio = currentHeight - heightN;
            if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                n.setParentDistance(p,n.getParentDistance(p) * ratio); // Coalescent_Unit = Tao /  theta/2
            }
        }

        return currentST;
    }

    // P(gt|ST)
    public List<Double> getGTSLLByST(List<UltrametricTree> geneTrees, Network currentST) {
        List<Double> gProST = new ArrayList<>();
        GeneTreeProbability geneTreeProbability = new GeneTreeProbability();

//        for(int i=0;i<geneTrees.size();i++){
//            Trees.autoLabelNodes((MutableTree)geneTrees.get(i).getTree());
//            gProST.add(getGTLLByST(geneTrees.get(i),currentST));
//        }
        gProST = geneTreeProbability.calculateGTDistribution(currentST,operator.getTreeByUTree(geneTrees),null,false);
        return gProST;
    }
    // P(gt|ST)
    public double getGTLLByST(UltrametricTree geneTree, Network<NetNodeInfo> currentST) {
        GeneTreeProbability geneTreeProbability = new GeneTreeProbability();

        GeneTreeBrSpeciesNetDistribution calculator = new GeneTreeBrSpeciesNetDistribution(currentST,null);
        return calculator.calculateGTDistribution((UltrametricTree) geneTree);
    }

    // P(gt|ST)
    public List<Double> getGTSLLBySTYF(List<UltrametricTree> geneTrees, Network<NetNodeInfo> currentST) {
        GeneTreeProbability geneTreeProbability = new GeneTreeProbability();
        List<Double> gProST = new ArrayList<>();

        double[] result = new double[geneTrees.size()];

        GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(currentST,operator.getTreeByUTree(geneTrees),null);

        g.calculateGTDistribution(result);
        //List<Tree> gts = new ArrayList<Tree>();

        for(int i = 0;i<result.length;i++)
            gProST.add(Math.log(result[i]));

        return gProST;
    }

    // P(Seq|gt) NOTICE: divergence time
    public double getSeqLLByGT(UltrametricTree gt, Alignment aln){
        Utils._SUBSTITUTION_MODEL = "JC";
        Frequencies freq = new Frequencies(aln, false);
        Trees.autoLabelNodes((MutableTree) gt.getTree());
        SubstitutionModel subst = new JukesCantor(freq);
        BeagleTreeLikelihood beagle = new BeagleTreeLikelihood(aln, gt, new SiteModel(subst), null);
        return beagle.calculateLogP();

    }

    public double getGTSLL(List<UltrametricTree> currentGeneSet, Network<NetNodeInfo> currentST) {
        double ll = 0.0;
        List<Double> gProST = getGTSLLBySTYF(currentGeneSet, currentST);
        for(int i = 0; i < currentGeneSet.size(); i++) {
            UltrametricTree gt = currentGeneSet.get(i);
            double p2 = getSeqLLByGT(gt, trueSeq.get(i));
            double p = (double) gProST.get(i) + p2;

            if (p2 > 0)
                System.out.println("p2: " + p2);

            ll += p;
        }
        return ll;
    }


    //IMPORTANT: ST : tap -> cUnit -> tao    ||  return GTS: tao
    public List<UltrametricTree> simGTSByMS(Network<NetNodeInfo> currentST, double halfTheta, int numGTs){
        //Scale ST
        Iterator<NetNode<NetNodeInfo>> it = currentST.bfs().iterator();
        //TODO scaleST
//        while(it.hasNext()) {
//            NetNode n = it.next();
//            if( n.getParents().iterator().hasNext()){
//                NetNode p = (NetNode) n.getParents().iterator().next();
//                if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
//                    n.setParentDistance(p, n.getParentDistance(p) / halfTheta); // Coalescent_Unit = Tao /  theta/2
//                }
//            }
//        }
        List<Tree> gts = simulator.simulateGeneTrees(currentST,null, refineGTSize);
        List<UltrametricTree> uGTS = operator.getUTreeByTree(gts);

        //Iterator uit = uGTS.iterator();

        //TODO scaleGT
//        while(uit.hasNext()) {
//            UltrametricTree thisGT = (UltrametricTree) uit.next();
//            Iterator itNode = thisGT.getTree().getNodes().iterator();
//            while (itNode.hasNext()) {
//                STINode n = (STINode) itNode.next();
//                //n.getParentDistance()
//                if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
//                    n.setParentDistance(n.getParentDistance() * halfTheta); // Coalescent_Unit = Tao /  theta/2
//                }
//            }
//        }


//        Iterator<NetNode<NetNodeInfo>> it1 = currentST.bfs().iterator();
//        while(it1.hasNext()) {
//            NetNode n = it1.next();
//            if( n.getParents().iterator().hasNext()){
//                NetNode p = (NetNode) n.getParents().iterator().next();
//                if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
//                    n.setParentDistance(p, n.getParentDistance(p) * halfTheta); // Coalescent_Unit = Tao /  theta/2
//                }
//            }
//        }
        return uGTS;
    }

//    public List<UltrametricTree> cUnitToTime(List<UltrametricTree> inputGTS){
//        Iterator it = inputGTS.iterator();
//        while(it.hasNext()){
//            UltrametricTree thisGT = (UltrametricTree) it.next().;
//            Iterator itNode = thisGT.getNodes().iterator();
//            while(itNode.hasNext()) {
//                NetNode n = (NetNode)itNode.next();
//                if( n.getParents().iterator().hasNext()){
//                    NetNode p = (NetNode) n.getParents().iterator().next();
//                    if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
//                        n.setParentDistance(p, n.getParentDistance(p) * halfTheta); // Coalescent_Unit = Tao /  theta/2
//                    }
//                }
//            }
//        }
//        return inputGTS;
//    }


    //scale all branches by random[0.9,1.0]
    public double scaleAll(List<UltrametricTree> gts, Network<NetNodeInfo> currentST,double preLL){
        //double preLL = getGTSLL(gts, currentST);
        double ratio = Math.random();
        //Tuning parameter
        ratio = ratio*0.05 + 0.975;

        operator.scaleSTNet(currentST, ratio,false);
        operator.scaleGTS(gts, ratio, false);
        gtsHLL[itNum] = Double.valueOf(currentHLL);
        //6.10 TODO: if full or half
        double aftLL = computeFullLL( currentST, gts);

        if(preLL>aftLL){
            operator.scaleSTNet(currentST, ratio,true);
            operator.scaleGTS(gts, ratio, true);
            System.out.println("------Not scale by : " + String.valueOf(ratio));
            gtsLL[itNum] = preLL;
            maxLL = preLL;
            finalGTS = gts;
            bestEstimatedST = currentST;
            return preLL;
        }

        else{
            System.out.println("Scale by " + String.valueOf(ratio));
            finalGTS = gts;
            maxLL = aftLL;
            gtsLL[itNum] = preLL;
            gtsHLL[itNum] = currentHLL;
            bestEstimatedST = currentST;
            return aftLL;
        }
    }

    //scale all branches by random[0.9,1.0]
    public double scaleAllHalf(List<UltrametricTree> gts, Network<NetNodeInfo> currentST,double preLL){
        double ratio = Math.random();
        ratio = ratio*0.05 + 0.975;

        operator.scaleSTNet(currentST, ratio,false);
        operator.scaleGTS(gts, ratio, false);
        gtsLL[itNum] = currentLL;
        System.out.println("gtsLL before compute() : " + gtsLL[itNum]);
        double aftLL = computeHalfLL( currentST, gts);
        System.out.println("gtsLL after compute(): " + currentLL);
        if(preLL>aftLL){
            operator.scaleSTNet(currentST, ratio,true);
            operator.scaleGTS(gts, ratio, true);
            System.out.println("------Not scale by : " + String.valueOf(ratio));
            gtsHLL[itNum] = preLL;
            finalGTS = gts;
            bestEstimatedST = currentST;
            return preLL;

        }

        else{
            System.out.println("Scale by " + String.valueOf(ratio));
            finalGTS = gts;
            gtsHLL[itNum] = aftLL;
            gtsLL[itNum] = computeFullLL( currentST, gts);
            maxLL = aftLL;
            return aftLL;
        }

        //gtsLL[itNum] = preLL;
        //maxLL = preLL;
        //finalGTS = gts;
        //bestEstimatedST = currentST;

    }

    public List<UltrametricTree> copyGTS(List<UltrametricTree> gts) {
        List<UltrametricTree> copy = new ArrayList<UltrametricTree>();
        for(UltrametricTree t : gts){
            copy.add(new UltrametricTree(t));
        }
        return copy;
    }
    public static void main(String[] args) throws IOException, ParseException, InterruptedException {
        IterativeImprovement ii = new IterativeImprovement();
        ii.iigt(trueSeq, ITERATION);
        InferOperator op = new InferOperator("/Users/doriswang/PhyloNet/Data/IIG/", "/Users/doriswang/PhyloNet/Data/IIG/output",0);
        op.load17Data("/Users/doriswang/PhyloNet/Data/17-taxon/004/ST0/1/", 1000, 32, trueSeq, new ArrayList<String>());

//        String stString = "(Marmoset,((Tarsier,((Galago,Mouse_Lemur)1:2.4018020540034675,((Horse,Sloth)1:2.6056351956741435,(Tree_Shrew,(Rat,Rabbit)1:0.8024369743576611)0.92:0.08540838439417373)1:1.9852014529353605)1:0.655284855661726)1:4.260329699696363,(Macaque,(Orangutan,((Human,Chimpanzee)1:0.6408490436022495,Gorilla)1:2.3251916476229617)1:2.5898315305010144)1:2.7005416959060815));";
//        //String pattern1 = "^\\d+";
//        //TODO: if outgroup cannot be recognised -> wrong
//        String pattern = "\\d+(\\.\\d+)?:";
//        Pattern p = Pattern.compile(pattern);
//        StringBuffer sb = new StringBuffer();
//        Matcher m = p.matcher(stString);
//        boolean hasStart = false;
//        int firstEnd = -1;
//        String temp = "";
//        int count = 1;
//        while (m.find()) {
//            m.appendReplacement(sb, "I" + count + ":");
//            count++;
//
//        }
//        m.appendTail(sb);
//
//        System.out.println(sb.toString());
//
//        Tree t = Trees.readTree(sb.toString());
//        Trees.autoLabelNodes((MutableTree) t);
//        Trees.rootAndRemoveOutgroup((STITree) t, "Marmoset");
//        System.out.println(t.toString());
//                    Trees.scaleBranchLengths((STITree) t, halfTheta);
//                    bw2.write(t.toString() + "\n");
    }


////
//        //BufferedReader out = new BufferedReader(new FileReader(outputFile));
////        String stString = "(A,((B,C)0.46:0.04652001563489282,((P,(N,O)0.68:0.20067069546215124)1:1.9924301646902058,(((D,E)1:1.9924301646902058,(F,(G,H)1:1.299282984130261)0.97:0.6061358035703156)0.98:0.7114963192281418,((L,M)1:1.9924301646902058,(K,(I,J)1:1.9924301646902058)0.97:0.6061358035703156)0.96:0.6061358035703156)0.97:0.6061358035703156)1:1.9924301646902058));";
////        //System.out.println(stString);
////        ArrayList<Double> paras = new ArrayList<Double>();
////        String[] s1 = stString.split(":");
////        String temp = "";
////        String topo = "";
////        String firstNode = "";
////        for(int i = 0; i<s1.length; i++){
////            if(i==s1.length-1){
////                int index1 = s1[i].indexOf(")");
////                paras.add(Double.valueOf(s1[i].substring(0,index1)));
////                topo += s1[i].substring(index1);
////                continue;
////            }
//////            if(i==0)
//////                firstNode = s1[i].split("'")[0];
////            int index1 = s1[i].indexOf(")");
////            if(index1==s1[i].length()-1) {
////                topo += s1[i];
////                continue;
////            }
////            if(index1 != -1) {
////                s1[i] = s1[i].substring(0, index1+1);
////                int index2 = s1[i].indexOf(",");
////                if (index2 == -1) {
////                    index2 = s1[i].indexOf(")");
////                    paras.add(Double.valueOf(s1[i].split("\\)")[0]));
////                    temp = s1[i].substring(index2);
////                    s1[i] = temp;
////                } else {
////                    if(i!=0) {
////                        paras.add(Double.valueOf(s1[i].substring(0, index2)));
////                        //temp = s1[i].split(",");
////                        temp = s1[i].substring(index2);
////                        s1[i] = temp;
////                    }
////                }
////            }
////            topo+=s1[i];
////        }
////
////        BniNetwork st = (BniNetwork) Networks.readNetwork(topo);
////        Networks.autoLabelNodes(st);
////        Networks.autoLabelNodes(st);
////        int paraLen = paras.size();
////        String stNet = st.toString();
////        String pattern = "I\\d+";
////        Pattern p = Pattern.compile(pattern);
////        Matcher m = p.matcher(stNet);
////        ArrayList<Integer> iNodes = new ArrayList<Integer>();
////        int count = 1;
////        while(m.find()){
////
////            String name = stNet.substring(m.start(),m.end());
////            if(name.equals("I1")||name.equals("I0"))
////                continue;
////            NetNode thisNode = st.findNode(name);
////            thisNode.setParentDistance((NetNode)thisNode.getParents().iterator().next(),paras.get(paraLen-count));
////            count++;
////        }
////        paras.get(3);
////        System.out.println(st.toString());
////        //TODO : just for YEAST
////        NetNode root = (NetNode) st.findNode("O").getParents().iterator().next();
////        st.resetRoot(root);
//        //out.close();
//    }
//

}
