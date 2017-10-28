package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

/**
 * Created by doriswang on 6/26/17.
 * For Snack data
 */


import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.likelihood.BeagleTreeLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.sitemodel.SiteModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.Frequencies;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.JukesCantor;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.TemporalConstraints;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GLASSInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.jcp.xml.dsig.internal.dom.DOMUtils;
import sun.nio.ch.Net;

import java.io.*;
//import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SYREstimator {
    protected static String _msdir;
    protected static String _ASTRALdir;
    private static String TREE_DIR;
    private static String SEQ_DIR;
    private static String RESULT_DIR;
    private static boolean TEST;

    private Network<NetNodeInfo> INIT_ST; //cUnit
    private Network<NetNodeInfo> currentST;
    private String GLOBALBESTST;
    private String trueST;
    //7.26
    private List<String> checkedTopo;  //topos as checked order
    private List<List<Tree>> optGTSList;  // For one locus l, locus<map<bestGT,LL>>


    private List<String> trueGTS;
    private List<Tree> trueTopos;
    private List<Tree> finalGTS;
    private static HashMap<String, Double> cHeights;
    private double[] ogHeight;

    private List<String> GLOBALBESTGTS;

    private static List<Alignment> trueSeq;
    private double GLOBALBESTLL;
    private IIGTSimulator simulator;
    private InferOperator operator;

    //TODO: PARA
    private static int ITERATION; //{1k, 5k ,10k}
    private static boolean ISREALDATA;
    private static boolean FULL_LL;
    private static String DATASET;
    private static double halfTheta;//need to init in ALL classes ||  theta -> POP_SIZE
    private static Integer lociNum;
    private static Integer refineGTSize;

    private static double[] _scales; //temp never used
    private static int[] _seqLens;  // 665,525,420,260,260,469,796,296,220,262,256,194,447,849,798,522,267,274,471,
    private static int taxaNum;
    private double[] gtsLL;// LL(Seq|GTS)*(GTS|ST)
    private double[] gtsHLL;
    private double[] currentITLL;
    private double currentLL;// LL(Seq|GTS)*(GTS|ST)
    private double currentHLL;

    private double[] currentGTLL;
    private double[] currentHGTLL;
    private double[] gtDis;
    private List<Double> bestSTLL;
    private List<Double> bestSTD;
    private List<Integer> bestITNum;
    private int bestLLIter;
    private double[][] gtDis_Locus;
    private double[][] msDist_Locus;

    private double[] stDis;
    //private double[][] stDis_Y;
    private double maxLL; //current best LL
    private double iter_LL;
    private double ll_Ratio;
    private int itNum;
    private long[] time;
    private static String STRATEGY; //RANDOM IMPROVE NONE
    private boolean ifOutGroup;

    //All parameters should be initialed here
    SYREstimator() throws IOException, InterruptedException, ParseException {

        //TEST = true;

        _msdir = "/Users/doriswang/PhyloNet/tools/msFiles/msdir/ms";
        _ASTRALdir = "/Users/doriswang/PhyloNet/tools/Astral/";
        TREE_DIR = "/Users/doriswang/PhyloNet/Data/IIG/tree/";
        SEQ_DIR = "/Users/doriswang/PhyloNet/Data/IIG/seq/";
        RESULT_DIR = "/Users/doriswang/PhyloNet/Data/IIG/result/";

        trueGTS = new ArrayList<String>();
        finalGTS = new ArrayList<Tree>();
        trueTopos = new ArrayList<Tree>();

        GLOBALBESTGTS = new ArrayList<String>();
        trueSeq = new ArrayList<Alignment>();
        //TODO: PARA
        ITERATION = 50; //{1k, 5k ,10k}
        ISREALDATA = true;
        FULL_LL = true;
        DATASET = "SYR";
        halfTheta = 0.001;//   1/2(10*e-2, 10*e-3)  0.0015(K) 0.002(Yang) need to init in ALL classes ||  theta -> POP_SIZE
        lociNum = 38;
        ifOutGroup = false;

        ogHeight = new double[2]; // og height from all initGT[min, max]
        refineGTSize = 30;
        _scales = new double[]{1.0};
        _seqLens = new int[]{665,525,420,260,260,469,796,296,220,262,256,194,447,849,798,522,267,274,471};
        taxaNum = 20;
        ll_Ratio = 1.0; // P(seq|GT): p(GT|ST)
        gtsLL = new double[ITERATION];// LL(Seq|GTS)*(GTS|ST)
        // LL(Seq|GTS)*(GTS|ST) for pseudo version || not used for classic version
        gtsHLL = new double[ITERATION];
        currentITLL = new double[ITERATION];
        currentLL = 0.0; // LL(Seq|GTS)*(GTS|ST)
        currentHLL = 0.0;
        iter_LL = 0.0;
        currentGTLL = new double[lociNum]; // P(Seq|GT)
        currentHGTLL = new double[lociNum];
        currentITLL = new double[ITERATION];
        cHeights = new HashMap<String,Double>();
        gtDis = new double[ITERATION];
        stDis = new double[ITERATION];
        //stDis_Y = new double[ITERATION][4];
        gtDis_Locus = new double[ITERATION][lociNum];
        msDist_Locus = new double[ITERATION][lociNum];
        time = new long[2];
        bestLLIter = 0;
        bestITNum = new ArrayList<Integer>();
        bestSTLL = new ArrayList<Double>();
        bestSTD = new ArrayList<Double>();
        //checkedTopo = new ArrayList<String>();  //topos as checked order
        //optGTSList = new ArrayList<List<Tree>>() ;
        STRATEGY = "RANDOM"; //RANDOM IMPROVE N
        //TODO Simulation
        simulator = new IIGTSimulator(lociNum, _scales, _seqLens, halfTheta, ITERATION);
        operator = new InferOperator(ITERATION);
        String resultFolder = RESULT_DIR + ITERATION + "/" + taxaNum + "_" + lociNum + "/" + refineGTSize + "/";
        operator.isExitsPath(resultFolder);
        if (ISREALDATA) {
            if (DATASET.equals(("SYR"))) {
                trueST = "((((((8,5),2),(4,1)),(7,(6,3))),((((10,11),9),(13,12)),(14,15))),((18,17),((19,20),16)));";
//                trueGTS.add("(1:0.03798,(((((2:0.004797,9:0.003204,3:0.000806,4:0.000772,5:0.000810,6:0.000787,7:0.000786,10:0.000798):0.001607,8:0.002159):0.037802,(((11:0.000817,13:0.000808,14:0.000809):0.003299,(12:0.000806,15:0.000804,16:0.000819,17:0.000806):0.002343):0.002431,(18:0.000788,19:0.000807):0.001933):0.029201):0.032477,((20:0.000816,(21:0.000801,22:0.000811,23:0.000801):0.004790,24:0.000806):0.001593,25:0.001260):0.049365):0.076642):0.031080);");
//                trueGTS.add("(1:0.03798,(((((2:0.004797,9:0.003204,3:0.000806,4:0.000772,5:0.000810,6:0.000787,7:0.000786,10:0.000798):0.001607,8:0.002159):0.037802,(((11:0.000817,13:0.000808,14:0.000809):0.003299,(12:0.000806,15:0.000804,16:0.000819,17:0.000806):0.002343):0.002431,(18:0.000788,19:0.000807):0.001933):0.029201):0.032477,((20:0.000816,(21:0.000801,22:0.000811,23:0.000801):0.004790,24:0.000806):0.001593,25:0.001260):0.049365):0.076642):0.031080);");
                List<Alignment> fullAln = operator.loadSYRByRX(lociNum,_seqLens,taxaNum+1,trueSeq);
                //TODO 9-11 load startTree?
                for(int i = 0;i<lociNum;i++){
                    trueGTS.add("((((((8,5),2),(4,1)),(7,(6,3))),((((10,11),9),(13,12)),(14,15))),((18,17),((19,20),16)));");
                }
//                for(int i = 0; i<lociNum; i++){
//                    BufferedWriter phy = new BufferedWriter(new FileWriter("/Users/doriswang/PhyloNet/tools/standard-RAxML-master/" + String.valueOf(i) + "/dna.phy"));
//                    Alignment a = fullAln.get(i);
//                    phy.write( "26 " + _seqLens[0] + '\n');
//                    phy.flush();
//                    List<String> names = a.getTaxaNames();
//                    Map<String, String> thisAln = a.getAlignment();
//                    for (int j = 0; j < 26; j++) {
//                        String name = names.get(j);
//                        String seq = thisAln.get(name);
//                        phy.write(name + '\t' + seq + '\n');
//                    }
//                    phy.flush();
//                    phy.close();
//                }
            }
        }
    }

    //  Input: alignments & t
    //  Output: GTs and STs
    public List<Tree> iigt(List<Alignment> aln, int t) throws IOException, ParseException, InterruptedException {
        //CHANGE: PARA
        itNum = 0;
        ITERATION = t;
        long start = System.currentTimeMillis();
        for(int i = 0; i<trueGTS.size();i++){
            trueTopos.add(Trees.readTree(operator.removeOutgroup(trueGTS.get(i))));
        }
        //Initialization
        Network<NetNodeInfo> tempST = initST(aln);
        //optimizeBL()
        itNum++;
        for (int i = 1; i < t; i++) {
            System.out.println("\n" + "#" + itNum + " iteration start! ");
            System.out.println(" Current ST is " + tempST.toString());
            refineGeneSet(tempST);
            itNum++;
        }
        long end = System.currentTimeMillis();
        long costtime = end - start;
        String resultFolder = RESULT_DIR + ITERATION + "/" + taxaNum + "_" + lociNum + "/" + refineGTSize + "/";
        operator.isExitsPath(resultFolder);
        BufferedWriter llOut1 = new BufferedWriter(new FileWriter(resultFolder + "RunningTime.txt"));
        llOut1.write("Running time:" + String.valueOf(costtime) + "\n");
        llOut1.write("Running time for p(GTS|ST):" + String.valueOf(time[0]) + "\n");
        llOut1.write("Running time for p(Seq|GTS)*p(GTS|ST):" + String.valueOf(time[1]) + "\n");

        // System.out.println("------Running time: " + costtime);
        operator.isExitsPath(resultFolder);
        llOut1.flush();
        llOut1.close();
        BufferedWriter llOut = new BufferedWriter(new FileWriter(resultFolder + "HLikelihood.txt"));
        for (int k = 0; k < ITERATION; k++) {
            llOut.write(k + " " + String.valueOf(gtsHLL[k]) + "\n");
        }
        //TODO llOut.write("AVE: " + operator.getAverage(gtsHLL));
        //TODO llOut.write("SD: " + operator.getStandardDevition(gtsHLL));
        llOut.flush();
        llOut.close();
        BufferedWriter llPOut = new BufferedWriter(new FileWriter(resultFolder + "FullLikelihood.txt"));
        gtsLL[0] = gtsLL[1];
        for (int k = 0; k < ITERATION; k++) {
            llPOut.write(k + " " + String.valueOf(gtsLL[k]) + "\n");
        }
        //TODO llPOut.write("AVE: " + operator.getAverage(gtsLL));
        //TODO llPOut.write("SD: " + operator.getStandardDevition(gtsLL));
        llPOut.flush();
        llPOut.close();
        for (int tid = 0; tid < lociNum; tid++) {

            String fileName = resultFolder + "GTDistance_" + tid + ".txt";
            BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
            for (int it = 0; it < ITERATION; it++) {
                out.write(it + " " + gtDis_Locus[it][tid] + "\n");
                out.flush();
            }
            out.close();
            String fileName1 = resultFolder + "MSGTDistance_" + tid + ".txt";
            BufferedWriter out1 = new BufferedWriter(new FileWriter(fileName1));
            for (int it = 0; it < ITERATION; it++) {
                out1.write(it + " " + msDist_Locus[it][tid] + "\n");
                out1.flush();
            }
            out1.close();
        }
        //String g_resultFolder = RESULT_DIR + ITERATION + "/global/";
        BufferedWriter bw2 = new BufferedWriter(new FileWriter(resultFolder + "GlobalBestST" + ".txt"));
        bw2.write(GLOBALBESTST + "\n");
        bw2.write("STATEGY : " + STRATEGY + "\n");

        bw2.write("Best Likelihood value : " + String.valueOf(GLOBALBESTLL) + "\n");
        bw2.write("Iteration number of best Likelihood value : " + String.valueOf(bestLLIter) + "\n");
        bw2.write("Likelihood combined ratio : " + ll_Ratio + "\n");

        bw2.write("Number of loci : " + lociNum + "\n");
        bw2.write("Number of refine gene tree set :  " + refineGTSize + "\n");
        bw2.write("Number of iteration :  " + ITERATION + "\n");
        bw2.write("SeqLens :  " + _seqLens[0] + "\n");
        bw2.write("Is real data :  " + ISREALDATA + "\n");
        bw2.write("Is full likelihood  " + FULL_LL + "\n");
        bw2.write("HalfTheta  " + halfTheta + "\n");

        bw2.flush();
        bw2.close();
        String fileName = resultFolder + "GlobalBestTrees.txt";
        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
        out.write(GLOBALBESTST);
        for (int tid = 0; tid < GLOBALBESTGTS.size(); tid++) {
            String gt = GLOBALBESTGTS.get(tid);
            out.write(gt + "\n");
            out.flush();
        }
        out.close();
        BufferedWriter ldOut = new BufferedWriter(new FileWriter(resultFolder + "distance.txt"));
        //6.10 TODO multi-ST distance
        double totalD = 0.0;
        System.out.println("---------------------------Result-----------------------------");
        System.out.println("Running time:" + String.valueOf(costtime));
        System.out.println("Best Likelihood value : " + String.valueOf(GLOBALBESTLL) );
        String stDistance = "RF of Best ST : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
        ldOut.write(stDistance + "\n");
        System.out.println(stDistance);
        String gtDistance = "RF of GT_";
        double thisD = 0.0;
        for (int k = 0; k < lociNum; k++) {
            thisD = operator.getDistance(Trees.readTree(GLOBALBESTGTS.get(k)), Trees.readTree(trueGTS.get(k)));
            ldOut.write(gtDistance + k + " :" + thisD + "\n");
            totalD += thisD;
        }
        ldOut.write("Total distance of GTS is " + totalD + "\n");
        ldOut.write("Average distance of GTS is " + totalD / trueGTS.size() + "\n");

        ldOut.write("---------------------------Global best ST change -----------------------------" + "\n");
        for(int i = 0; i< bestSTD.size();i++){
            ldOut.write(String.valueOf(bestITNum.get(i) + " : " + String.valueOf(bestSTD.get(i) + "  :  " + String.valueOf(bestSTLL.get(i))))+ "\n");
        }

        System.out.println("Average distance of GTS is " + totalD / trueGTS.size());
        System.out.println("Best ST is " + GLOBALBESTST);
        BufferedWriter gtD = new BufferedWriter(new FileWriter(resultFolder + "GTdistance.txt"));
        for (int gtn = 0; gtn < ITERATION; gtn++) {
            gtD.write(String.valueOf(gtn) + ":" + String.valueOf(gtDis[gtn]) + "\n");
            //TODO 8-10
        }
        gtD.flush();
        gtD.close();
        ldOut.flush();
        ldOut.close();


        if (DATASET.equals("YEAST")) {
//            for (int i = 0; i < 4; i++) {
//                BufferedWriter stD = new BufferedWriter(new FileWriter(resultFolder + "STdistance_" + String.valueOf(i) + ".txt"));
//                for (int stn = 0; stn < ITERATION; stn++) {
//                    stD.write(String.valueOf(stn) + ":" + String.valueOf(stDis_Y[stn][i]) + "\n");
//                }
//                stD.flush();
//                stD.close();
//            }
        } else {

            BufferedWriter stD = new BufferedWriter(new FileWriter(resultFolder + "STdistance.txt"));
            for (int stn = 0; stn < ITERATION; stn++) {
                stD.write(String.valueOf(stn) + ":" + String.valueOf(stDis[stn]) + "\n");
            }
            stD.flush();
            stD.close();
//            BufferedWriter sFtLL = new BufferedWriter(new FileWriter(resultFolder + "STITLL.txt"));
//            for (int stn = 0; stn < ITERATION; stn++) {
//                sFtLL.write(String.valueOf(stn) + " " + String.valueOf(iter_LL) + "\n");
//            }
//
//            sFtLL.flush();
//            sFtLL.close();

            BufferedWriter itLL = new BufferedWriter(new FileWriter(resultFolder + "STCurrentLL.txt"));
            for (int stn = 0; stn < ITERATION; stn++) {
                itLL.write(String.valueOf(stn) + " " + String.valueOf(currentITLL[stn]) + "\n");
            }

            itLL.flush();
            itLL.close();

            llOut = new BufferedWriter(new FileWriter(resultFolder + "LL_Analysis.txt"));
            for (int k = 0; k < ITERATION; k++) {
                llOut.write(k + " " + String.valueOf(gtsHLL[k]) + "\n");
            }
            llOut.write("Half AVE: " + operator.getAverage(gtsHLL) + "\n");
            llOut.write("Half SD: " + operator.getStandardDevition(gtsHLL) + "\n");

            llOut.write("Full AVE: " + operator.getAverage(gtsLL) + "\n");
            llOut.write("Full SD: " + operator.getStandardDevition(gtsLL) + "\n");
            llOut.write("Current AVE: " + operator.getAverage(currentITLL) + "\n");
            llOut.write("Current SD: " + operator.getStandardDevition(currentITLL) + "\n");
            llOut.flush();
            llOut.close();
        }

        return finalGTS;
    }

    //  Input:  tempBestST
    //  Output: current best GTS
    public List<Tree> refineGeneSet (Network < NetNodeInfo > tempST) throws ParseException, IOException {
        List<Tree> currentTopo = simulator.generateGTS(tempST,refineGTSize, _msdir,null, cHeights, false); //operator.MS
        operator.getMSDist(trueTopos, currentTopo,msDist_Locus, itNum);
        List<Tree> ogTopos = new ArrayList<Tree>();
        for(int i = 0; i< currentTopo.size(); i++){
            String tStr = currentTopo.get(i).toString();
            String newT = "(" + tStr.substring(0,tStr.length()-1) + ",O);";
            MutableTree tempT = (MutableTree)Trees.readTree(newT);
            Trees.autoLabelNodes(tempT);
            ogTopos.add(tempT);
        }

        List<Tree> gts = getBestGTSByR(ogTopos);
        boolean move = false;
        //double p2 = currentLL - currentHLL;
        currentITLL[itNum] = currentLL;
        iter_LL = currentLL;
        List<Tree> gtsOG = new ArrayList<Tree>();
        for(int i = 0; i<gts.size(); i++){
            gtsOG.add(operator.rerootRAxML(gts.get(i),"O"));
        }
        gts = operator.getNoOGGTS(gtsOG);
        Tree inferredST = operator.inferSTByGLASS(gtsOG);
        if(inferredST.getNode("O").getParentDistance()!=Double.NEGATIVE_INFINITY)
            gts.add(inferredST);
//        if (FULL_LL)
//            iter_LL = currentHLL*ll_Ratio + p2*(2-ll_Ratio);
//        else
//            iter_LL = currentHLL;

        //currentLL ->  iter_LL
        if (STRATEGY == "NONE") {
            move = true;
        } else if (STRATEGY == "IMPROVE") {
            if (iter_LL > maxLL) {
                move = true;
                maxLL = iter_LL;
                bestLLIter = itNum;
                GLOBALBESTST = tempST.toString();
                bestSTLL.add(iter_LL);
                bestSTD.add(operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)));
                bestITNum.add(itNum);
                System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GLOBALBEST ST change");
                System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ distance : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)));
                GLOBALBESTLL = iter_LL;
                Iterator it = gts.iterator();
                GLOBALBESTGTS = new ArrayList<String>();
                while (it.hasNext()) {
                    GLOBALBESTGTS.add(it.next().toString());

                }
            } else
                move = false;
        } else if (STRATEGY == "RANDOM") {
            if (iter_LL > maxLL) {
                move = true;
                if (iter_LL > GLOBALBESTLL) {
                    maxLL = iter_LL;
                    bestLLIter = itNum;
                    GLOBALBESTLL = iter_LL;
                    GLOBALBESTST = tempST.toString();
                    Iterator it = gts.iterator();
                    GLOBALBESTGTS = new ArrayList<String>();
                    while (it.hasNext()) {
                        GLOBALBESTGTS.add(it.next().toString());
                    }
                    System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GLOBALBEST ST change");
                    System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ distance : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)));
                }
            } else {
                double random = Math.random();
                if (random < 0.98)
                    move = false;
                else {
                    System.out.println("--------Random-Accept: #Iteration " + itNum + " : " + currentLL);
                    move = true;
                }
            }
        }
        if(itNum==1)
            move = true;

        if (move == true) {
//            List<Tree> gtsOG = new ArrayList<Tree>();
//            for(int i = 0; i<gts.size(); i++){
//                gtsOG.add(operator.rerootRAxML(gts.get(i),"O"));
//            }
//            gts = operator.getNoOGGTS(gtsOG);
//            Tree inferredST = operator.inferSTByGLASS(gts);
//            if(inferredST.getNode("O").getParentDistance()!=Double.NEGATIVE_INFINITY)
//                gts.add(inferredST);
            tempST = inferSTByAST(gtsOG);

            maxLL = iter_LL;
            //currentITLL[itNum] = iter_LL;
            stDis[itNum] = operator.getDistance(Trees.readTree(tempST.toString()), Trees.readTree(trueST));
            System.out.println("ST Distance of # " + itNum + " is : " + stDis[itNum]);

            double thisD = 0.0;
            for (int k = 0; k < lociNum; k++) {

                double temp = operator.getDistance(gts.get(k), Trees.readTree(trueGTS.get(k)));
                gtDis_Locus[itNum][k] = temp;
                thisD += temp;
                double tt = msDist_Locus[itNum][k];
            }
            gtDis[itNum] = thisD / lociNum;
            System.out.println("GTS Distance of # " + itNum + " is : " + gtDis[itNum]);
            System.out.println("**********Accept ST  : " + tempST.toString());
            gtsLL[itNum] = currentLL;
            gtsHLL[itNum] = currentHLL;

        } else {
            System.out.println("----------Reject ST  : " + tempST.toString());
            gtsLL[itNum] = gtsLL[itNum - 1];
            gtsHLL[itNum] = gtsHLL[itNum - 1];
            //currentITLL[itNum] = currentITLL[itNum-1];
            //6.11 TODO: scale by ll or half_LL
            stDis[itNum] = stDis[itNum - 1];
            System.out.println("ST Distance of # " + itNum + " is : " + stDis[itNum]);
            gtDis[itNum] = gtDis[itNum - 1];
            System.out.println("GTS Distance of # " + itNum + " is : " + gtDis[itNum]);

            for (int k = 0; k < lociNum; k++) {
                gtDis_Locus[itNum][k] = gtDis_Locus[itNum-1][k];
            }
        }
        System.out.println("Current combined ll : " + iter_LL);
        System.out.println("Full ll : " + gtsLL[itNum] + " ||| " + "Half ll : " + gtsHLL[itNum]+ " ||| " + "Combined ll : " + currentITLL[itNum]);
        return finalGTS;
    }


    //Input: topos, ST
    //Output: best p(seq|G) RAxML trees(17)
    public List<Tree> getBestGTSByR(List<Tree> topos) throws IOException {
        //double[][] ll = new double[lociNum][refineGTSize];
        currentHLL = 0.0;
        currentLL = 0.0;
        List<Tree> bestGTs = new ArrayList<Tree>();
        List<Double> bestGTLL = new ArrayList<Double>();
        List<Integer> bestGTNum = new ArrayList<Integer>();
        //List<List<Tree>> treeList  = new ArrayList<List<Tree>>();
        operator.updateMSTopos(topos);
        operator.runUpdateShell(refineGTSize,lociNum);
        for(int i = 0 ; i<lociNum; i++){
            double[] llList = operator.getLocusLL(i, refineGTSize);
            double bestLL = llList[0];
            int bestGTNumi = 0;
            for(int j = 0; j < llList.length; j++) {
                if (bestLL < llList[j]) {
                    bestLL = llList[j];
                    bestGTNumi = j;
                }
            }
            if(i==0||i==18)
                bestGTs.add(operator.scaleGT(operator.getbestGTi(i,bestGTNumi),halfTheta*0.5,false)); // mtDNA_theta = 0.5* general theta
            else
                bestGTs.add(operator.scaleGT(operator.getbestGTi(i,bestGTNumi),halfTheta,false));
            bestGTLL.add(bestLL);
            currentHLL += bestLL;
            bestGTNum.add(bestGTNumi);
        }
        currentLL = currentHLL;
        return bestGTs;
    }

    //Run updateShell
    //Input: refineSize, loci,
    //TODO 8 888 move to RMS
//    public List<Tree> optimizeBLGT(List<Tree> topos, int lociNum, int msSize, Double halfTheta) throws IOException {
//        updateMSTopos(topos);
//        runUpdateShell(msSize, lociNum);
//        List<Tree> bestGT = new ArrayList<Tree>();
//        for(int i = 0 ; i<lociNum; i++){
//            double[] llList = getLocusLL(i, msSize);
//            double bestLL = llList[0];
//            int bestGTNum = 0;
//            for(int j = 0; j < llList.length; j++) {
//                if (bestLL > llList[j]) {
//                    bestLL = llList[j];
//                    bestGTNum = j;
//                }
//            }
//            String gtFile = _RAxMLdir + "RAxML_result.L" + String.valueOf(i) + "T" + String.valueOf(bestGTNum);
//            BufferedReader gtReader = new BufferedReader(new FileReader(gtFile));
//            String gtString = gtReader.readLine().trim();
//            bestGT.add(Trees.readTree(gtString));
//            gtReader.close();
//        }
//        return  bestGT;
//    }

//
//    //Input: topos, ST
//    //Output: best p(seq|G) RAxML trees
//    public List<Tree> getBestGTSByR(List<Tree> topos, List<Alignment> alns) throws IOException {
//        double[][] ll = new double[lociNum][refineGTSize];
//        currentHLL = 0.0;
//        currentLL = 0.0;
//        List<Tree> bestGTs = new ArrayList<Tree>();
//        List<List<Tree>> treeList  = new ArrayList<List<Tree>>();
//        operator.updateMSTopos(topos);
//        for(int i = 0; i< lociNum; i++){
//            List<Tree> rTrees = operator.optimizeBLGT(topos,alns, halfTheta);
//            treeList.add(rTrees);
//            for(int j = 0; j<lociNum; j++){
//                Tree temp = rTrees.get(j);
//                ll[j][i] = temp.getRoot().getParentDistance();
//                temp.getRoot().setParentDistance(Double.NEGATIVE_INFINITY);
//            }
//        }
//        for(int k = 0; k<lociNum; k++){
//            double bestLL = Double.NEGATIVE_INFINITY;
//            Tree bestT = (Tree)treeList.get(0).get(0);
//            for(int i =0; i< refineGTSize; i++) {
//                if (ll[k][i] > bestLL) {
//                    bestLL = ll[k][i];
//                    bestT = treeList.get(i).get(k);
//                }
//            }
//            currentHGTLL[k] = bestLL;
//            currentGTLL[k] = currentHGTLL[k];
//            currentLL += bestLL;
//            currentHLL = currentLL;
//            bestGTs.add(bestT);
//        }
//
//        return bestGTs;
//    }

//    public List<Tree>  getSeqLL(List<Tree> topos, List<Alignment> alns) throws IOException{
//        double[][] ll = new double[lociNum][refineGTSize];
//        List<Tree> bestGTs = new ArrayList<Tree>();
//        List<List<Tree>> treeList  = new ArrayList<List<Tree>>();
//        for(int i = 0; i< refineGTSize; i++){
//            List<Tree> rTrees = operator.optimizeBLGT(topos.get(i).toNewick(),alns, halfTheta);
//            treeList.add(rTrees);
//            for(int j = 0; j<lociNum; j++){
//                Tree temp = rTrees.get(j);
//                ll[j][i] = temp.getRoot().getParentDistance();
//                temp.getRoot().setParentDistance(Double.NEGATIVE_INFINITY);
//            }
//        }
//        for(int k = 0; k<lociNum; k++){
//            double bestLL = Double.NEGATIVE_INFINITY;
//            Tree bestT = (Tree)treeList.get(0).get(0);
//            for(int i =0; i< refineGTSize; i++) {
//                if (ll[k][i] > bestLL) {
//                    bestLL = ll[k][i];
//                    bestT = treeList.get(i).get(k);
//                }
//            }
//            currentHGTLL[k] = bestLL;
//                bestGTs.add(bestT);
//        }
//
//        return bestGTs;
//    }

    //Input: #locus, topos, aln_l
    //Output: currentGTll[i] = best p(seq|G) || RAxML tree for locus L
//    public Tree getBestGTByR(int locus, List<Tree> topos) throws IOException {
//        double bestLL = Double.NEGATIVE_INFINITY;
//        double bestHLL = Double.NEGATIVE_INFINITY;
//
//        Alignment aln = trueSeq.get(locus);
//        List<Tree> rTrees = operator.optimizeBLGT()
////        for (int i = 0; i < currentGeneSet.size(); i++) {
////            //List<Tree> tList1 = new ArrayList<Tree>();
////            //tList1.add(currentGeneSet.get(i));
////            Tree gt = operator.scaleGT(currentGeneSet.get(i), halfTheta, true);
////
////            double p2 = getSeqLLByGT(gt, aln);
////            if (p2 > 0) {
////                p2 = getSeqLLByGT(gt, aln);
////            }
////
////            if (p2 > 0) {
////                System.out.println("p(Seq|GT):" + p2 + "  || Locus: " + locus);
////                continue;
////            }
////            double p = gProST.get(i) + p2;
////            if (FULL_LL) {
////                if (p >= bestLL) {
////                    bestT = gt;
////                    bestLL = p;
////                    currentGTLL[locus] = bestLL;
////                    currentHGTLL[locus] = p2;
////
////                }
////            } else {
////                if (p2 >= bestLL) {
////                    bestT = gt;
////                    bestLL = p2;
////                    currentHGTLL[locus] = bestLL;
////                    currentGTLL[locus] = p;
////                }
////            }
////            operator.scaleGT(currentGeneSet.get(i), halfTheta, false);
////        }
//
//        return bestT;
//    }

    public List<Tree> getBestGTSByMS(List<Tree> currentGTS, Network tempST) throws IOException {
        List<Double> gProST = new ArrayList<>();
        gProST = getGTSLLBySTYF(currentGTS, tempST);
        List<Tree> gts = new ArrayList<>();
        double ll = 0.0;
        double hLL = 0.0;
        for (int i = 0; i < lociNum; i++) {
            Tree bestGTi = getBestGTByMS(i, currentGTS, gProST);
            gts.add(bestGTi);
            ll += currentGTLL[i];
            hLL += currentHGTLL[i];
        }
        currentLL = ll;
        currentHLL = hLL;
        return gts;
    }

    public Tree getBestGTByMS(int locus, List<Tree> currentGeneSet, List<Double> gProST) throws IOException {
        double bestLL = Double.NEGATIVE_INFINITY;
        double bestHLL = Double.NEGATIVE_INFINITY;

        Alignment aln = trueSeq.get(locus);
        Tree bestT = currentGeneSet.get(0);
        for (int i = 0; i < currentGeneSet.size(); i++) {
            //List<Tree> tList1 = new ArrayList<Tree>();
            //tList1.add(currentGeneSet.get(i));
            Tree gt = operator.scaleGT(currentGeneSet.get(i), halfTheta, true);

            double p2 = getSeqLLByGT(gt, aln);
            if (p2 > 0) {
                p2 = getSeqLLByGT(gt, aln);
            }

            if (p2 > 0) {
                System.out.println("p(Seq|GT):" + p2 + "  || Locus: " + locus);
                continue;
            }
            double p = gProST.get(i) + p2;
            if (FULL_LL) {
                if (p >= bestLL) {
                    bestT = gt;
                    bestLL = p;
                    currentGTLL[locus] = bestLL;
                    currentHGTLL[locus] = p2;

                }
            } else {
                if (p2 >= bestLL) {
                    bestT = gt;
                    bestLL = p2;
                    currentHGTLL[locus] = bestLL;
                    currentGTLL[locus] = p;
                }
            }
            operator.scaleGT(currentGeneSet.get(i), halfTheta, false);
        }

        return bestT;
    }


    // P(gt|ST)
    public List<Double> getGTSLLBySTYF(List<Tree> geneTrees, Network<NetNodeInfo> currentST) {
        List<Double> gProST = new ArrayList<>();
        double[] result = new double[geneTrees.size()];
        Network<NetNodeInfo> ultraST = operator.getUltraST(currentST.toString());
        GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(ultraST, geneTrees, null);
        g.calculateGTDistribution(result);
        for (int i = 0; i < result.length; i++)
            gProST.add(Math.log(result[i]));

        return gProST;
    }


    // Output: GTS with BL: mutations/site
    public Network<NetNodeInfo> initST(List<Alignment> aln) throws IOException, ParseException, InterruptedException {

        //List<Tree> tempOGTs = operator.initByRAxML(aln, _seqLens[0],taxaNum+1, lociNum);
        String[] strs = operator.getRAxMLStartT(lociNum); // strs[i] = treeString + " " + (Double) ll
        operator.getUpdateSh(lociNum,refineGTSize);
        String trees[]  = new String[lociNum];

//        trees[0] =    "(O:0.08930026338479349346,((((((14:0.00000100000050002909,13:0.00000100000050002909):0.00000100000050002909,11:0.00000100000050002909):0.00463757866709454004,(12:0.00000100000050002909,((17:0.00000100000050002909,16:0.00000100000050002909):0.00000100000050002909,15:0.00000100000050002909):0.00000100000050002909):0.00305993868333025512):0.00254876514318441065,(18:0.00000100000050002909,19:0.00000100000050002909):0.00206851264895196616):0.05651577787217316323,(8:0.00000100000050002909,((2:0.00151800868215477445,9:0.00000100000050002909):0.00304375182624630654,((5:0.00000100000050002909,((7:0.00000100000050002909,3:0.00000100000050002909):0.00000100000050002909,4:0.00000100000050002909):0.00000100000050002909):0.00000100000050002909,(10:0.00000100000050002909,6:0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00151677649814329686):0.10491674380054012705):0.08284202118961710160,(25:0.00000100000050002909,((24:0.00000100000050002909,20:0.00000100000050002909):0.00000100000050002909,(21:0.00000100000050002909,(23:0.00000100000050002909,22:0.00000100000050002909):0.00000100000050002909):0.00775837951559163550):0.00153472285054641008):0.16005822988799844886):0.33123363233198310773,1:0.07110686268781551966):0.0;";
//
//        trees[1] =     "(O:0.08930026338479349346,((((((14:0.00000100000050002909,13:0.00000100000050002909):0.00000100000050002909,11:0.00000100000050002909):0.00463757866709454004,(12:0.00000100000050002909,((17:0.00000100000050002909,16:0.00000100000050002909):0.00000100000050002909,15:0.00000100000050002909):0.00000100000050002909):0.00305993868333025512):0.00254876514318441065,(18:0.00000100000050002909,19:0.00000100000050002909):0.00206851264895196616):0.05651577787217316323,(8:0.00000100000050002909,((2:0.00151800868215477445,9:0.00000100000050002909):0.00304375182624630654,((5:0.00000100000050002909,((7:0.00000100000050002909,3:0.00000100000050002909):0.00000100000050002909,4:0.00000100000050002909):0.00000100000050002909):0.00000100000050002909,(10:0.00000100000050002909,6:0.00000100000050002909):0.00000100000050002909):0.00000100000050002909):0.00151677649814329686):0.10491674380054012705):0.08284202118961710160,(25:0.00000100000050002909,((24:0.00000100000050002909,20:0.00000100000050002909):0.00000100000050002909,(21:0.00000100000050002909,(23:0.00000100000050002909,22:0.00000100000050002909):0.00000100000050002909):0.00775837951559163550):0.00153472285054641008):0.16005822988799844886):0.33123363233198310773,1:0.07110686268781551966):0.0;";
        List tempGTs = new ArrayList<Tree>();
        double ogMin = Double.POSITIVE_INFINITY;
        double ogMax = 0.0;
        for(int i=0;i<lociNum;i++){
            trees = strs[i].split(" ");
            Tree tempT = Trees.readTree(trees[0]);
            tempT = operator.rerootRAxML(tempT,"O");
            //Tree tempT = tempOGTs.get(i);
            if(i==0||i==19)
                operator.scaleGT(tempT,halfTheta*0.5,false); // mtDNA_theta = 0.5* general theta
            else
                operator.scaleGT(tempT,halfTheta,false);
            tempGTs.add(tempT);
            String temp = operator.rerootAndRemove(tempT.toString(),"O");
            finalGTS.add(Trees.readTree(temp));

            GLOBALBESTGTS.add(temp);
            double ogH = tempT.getNode("O").getParentDistance();
            if(ogH>ogMax)
                ogMax = ogH;
            if(ogH<ogMin)
                ogMin = ogH;
        }
        ogHeight[0] = ogMin;
        ogHeight[1] = ogMax;


//        currentHGTLL[0] = -2248.794461;
//        currentHGTLL[1] = -2248.794461;
        //List<Tree> gtsOG = operator.getOGGTS(gts, ogHeight[0]);
        INIT_ST = inferSTByAST(tempGTs);
        //optimizeBL(tempGTs,INIT_ST,ogHeight[0]);
        //List<Double> proGT = getGTSLLBySTYF(finalGTS, INIT_ST);
        currentHLL = 0.0;
        currentLL = 0.0;

        for(int i = 0;i<currentHGTLL.length;i++) {
            currentHLL += currentHGTLL[i];
            currentGTLL[i] = currentHGTLL[i];

        }
        currentLL = currentHLL;
        GLOBALBESTLL = currentHLL;

//        double p2 = currentHLL-currentLL;
//        if (FULL_LL)
//            iter_LL = currentHLL*ll_Ratio + p2*(2-ll_Ratio);
//        else
        iter_LL = currentHLL;

        gtsHLL[0] = currentHLL;
        gtsLL[0] = currentLL;
        GLOBALBESTST = INIT_ST.toString();
        //trueST = operator.rerootAndRemove(trueST,"O");
        stDis[0] = operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
        double thisD = 0.0;
        //finalGTS = tempGTs;
        for (int k = 0; k < lociNum; k++){
            //trueGTS.set(k,operator.rerootAndRemove(trueGTS.get(k),"O"));
            double d = operator.getDistance((Tree) finalGTS.get(k), Trees.readTree(trueGTS.get(k)));
            gtDis_Locus[0][k] = d;
            msDist_Locus[0][k] = 0;
            thisD += d;
        }
        gtDis[0] = thisD / lociNum;
        return INIT_ST;
    }


    public static void main(String[] args) throws IOException, ParseException, InterruptedException {

        SYREstimator syr = new SYREstimator();
        InferOperator operator = new InferOperator(ITERATION);
        IIGTSimulator simulator1 = new IIGTSimulator(2, _scales, _seqLens, halfTheta, ITERATION);
        //ast.initST(trueSeq);
        syr.iigt(trueSeq,ITERATION);
        System.out.println("Finish");
    }

    //input: rooted gts with OG
    //output: rooted ST without OG
    public BniNetwork<NetNodeInfo> inferSTByAST(List<Tree> gtsOG) throws IOException, ParseException {
        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i /Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";
        String cmdFile = _ASTRALdir + "tempCMD.sh";
        BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
        cmd.write(command + '\n');
        cmd.flush();
        cmd.close();
        String inputFile = _ASTRALdir + "test_data/tempIn.tre";
        BufferedWriter in = new BufferedWriter(new FileWriter(inputFile));
        Tree temp = gtsOG.get(0);

        for (int i = 0; i < gtsOG.size(); i++) {
            in.write(gtsOG.get(i).toString() + "\n");
        }
        in.flush();
        in.close();
        ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
        pb.redirectErrorStream(true);
        try {
            Process proc = pb.start();
            try {
                proc.waitFor();
            } catch (InterruptedException ex) {
                System.out.println(ex.getMessage());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        String outputFile = _ASTRALdir + "test_data/testOut.tre";
        BufferedReader out = new BufferedReader(new FileReader(outputFile));
        String stString = out.readLine().trim();
        int index0 = 0;
        int index1 = 1;
        List<Integer> indexes = new ArrayList<Integer>();
        for(int i = 1;i<stString.length();i++) {
            if (stString.charAt(index1) == ':') {
                indexes.add(i);
            }
            index1++;
        }
        for(int i = 1;i<=indexes.size();i++) {
            index1 = indexes.get(indexes.size() - i);
            for (int j = 1; j < stString.length(); j++) {
                index0 = index1 - j;
                if (stString.charAt(index0) == ')') {
                    break;
                } else
                    index0--;
            }
            String newST = stString.substring(0, index0 + 1) + "I" + Integer.toString(i - 1) + stString.substring(index1);
            stString = newST;
        }
        System.out.println("AST String :"+stString);

//TODO 8.3
        Network n = operator.reRootAST(Trees.readTree(stString) ,"O");
        optimizeBL(gtsOG,n,ogHeight[0]);
        return (BniNetwork<NetNodeInfo>) n;
    }


    //  for all branches in st, if length = 0 -> length = ogHeight/10000
    //  if gts are ultrametric, add glassTree to bound ST
    public Network<NetNodeInfo> optimizeBL(List<Tree> gtsOG, Network st, Double ogHeight) throws IOException, ParseException {

        Map<Integer, List<String>> map = (Map<Integer, List<String>>)operator.getLevels(st);
        Set keys = map.keySet();
        int[] key = new int[keys.size()];
        Iterator keyIt = keys.iterator();
        for(int i = 0; i< keys.size();i++){
            key[i] = (int)keyIt.next();
        }
        int temp = 0;
        //bubble sort
        for(int i = 0; i<key.length-1;i++){
            for(int j = 0;j<key.length-1;j++){
                if(key[j]<key[j+1]){
                    temp = key[j];
                    key[j] = key[j+1];
                    key[j+1] = temp;
                }
            }
        }

        HashMap<String,Double> height = new HashMap<String,Double>();
        for(int i =0; i<keys.size(); i++){
            List<String> nodes = map.get(key[i]);

            Iterator<String> nodeIt = nodes.iterator();
            while(nodeIt.hasNext()){
                BniNetNode tempNode = (BniNetNode)st.findNode(nodeIt.next());
                BniNetNode p = (BniNetNode) tempNode.getParents().iterator().next();
                if(tempNode.getParentDistance(p)!=Double.NEGATIVE_INFINITY)
                    continue;
                    //never checked
                else {
                    Iterator cIt = p.getChildren().iterator();
                    BniNetNode sib = (BniNetNode) cIt.next();
                    if (sib.getName() == tempNode.getName())
                        sib = (BniNetNode) cIt.next();
                    // sib is leaf
                    if (sib.isLeaf() == true) {
                        double distance = Double.POSITIVE_INFINITY;
                        distance = getDist(tempNode.getName(), sib.getName(), gtsOG, st);
                        if(distance==0||distance == Double.POSITIVE_INFINITY) {
                            distance = ogHeight / 5000;
                        }
                        tempNode.setParentDistance(p, distance / 2);
                        sib.setParentDistance(p, distance / 2);
                        height.put(p.getName(),distance/2);

                    }
                    //sib is non-leaf
                    else {
                        // lowerDist
                        if(sib.getParentDistance(p)==0||sib.getParentDistance(p)==Double.NEGATIVE_INFINITY)
                            sib.setParentDistance(p,ogHeight/10000);
                        if (height.containsKey(sib.getName())) {
                            tempNode.setParentDistance(p, height.get(sib.getName()) + sib.getParentDistance(p));
                            height.put(p.getName(), height.get(sib.getName()) + sib.getParentDistance(p));
                        } else {
                            Double sibHeight = operator.getInternalHeight(sib, st, height);
                            height.put(sib.getName(), sibHeight);
                            //getDist(((BniNetNode) children.next()).getName(), ((BniNetNode) children.next()).getName(), gts, st);
                            double distance = sibHeight + sib.getParentDistance(p);
                            tempNode.setParentDistance(p, distance);
                            height.put(p.getName(), distance);
                        }
                    }
                }
            }
        }

        return st;
    }

    //input: Leaf_1 Leaf_2 gts st  (Two node must be LEAF NODE!!)
    //output: the shortest distance in all gts.
    public double getDist(String nodeName1, String nodeName2, List<Tree> gts, Network st) {
        double distance = Double.POSITIVE_INFINITY;

        // getDist(nonChecked leafnode1, leafnode 2 ...)
        for (int i = 0; i < gts.size(); i++) {
            Tree gt = gts.get(i);
            Trees.autoLabelNodes((MutableTree) gt);
            TNode n1 = gt.getNode(nodeName1);
            TNode n2 = gt.getNode(nodeName2);
            TNode root = gt.getRoot();
            TNode lcp = operator.getLCP(root, n1, n2);
            if (lcp.isLeaf())
                System.err.println("Wrong LCP!!");
            //System.out.println(n1.getName()+" + " + n2.getName()+" -> "  + lcp.getName());
            double temp = operator.getLocalDist(n1, n2, lcp, gt);
            if (temp < distance)
                distance = temp;
        }
        if(distance==Double.POSITIVE_INFINITY)
            System.err.println("Wrong Distance!!" + nodeName1 + " + " + nodeName2);

        return distance;
    }

    public Tree getBestGT(int locus, List<Tree> currentGeneSet, List<Double> gProST) throws IOException {
        double bestLL = Double.NEGATIVE_INFINITY;
        double bestHLL = Double.NEGATIVE_INFINITY;

        Alignment aln = trueSeq.get(locus);
        Tree bestGeneTree = currentGeneSet.get(0);
        for (int i = 0; i < currentGeneSet.size(); i++) {
            //List<Tree> tList1 = new ArrayList<Tree>();
            //tList1.add(currentGeneSet.get(i));
            Tree gt = operator.scaleGT(currentGeneSet.get(i), halfTheta, false);

            double p2 = getSeqLLByGT(gt, aln);
            if (p2 > 0) {
                p2 = getSeqLLByGT(gt, aln);
            }

            if (p2 > 0) {
                System.out.println("p(Seq|GT):" + p2 + "  || Locus: " + locus);
                continue;
            }
            double p = gProST.get(i) + p2;
            if (FULL_LL) {
                if (p >= bestLL) {
                    bestGeneTree = gt;
                    bestLL = p;
                    currentGTLL[locus] = bestLL;
                    currentHGTLL[locus] = p2;

                }
            } else {
                if (p2 >= bestLL) {
                    bestGeneTree = gt;
                    bestLL = p2;
                    currentHGTLL[locus] = bestLL;
                    currentGTLL[locus] = p;
                }
            }
            operator.scaleGT(currentGeneSet.get(i), halfTheta, true);
        }

        return bestGeneTree;
    }

    // P(gt|ST)
    public List<Double> getGTSLLByST(List<Tree> geneTrees, Network currentST) {
        GeneTreeProbability geneTreeProbability = new GeneTreeProbability();
        //TODO:if YF
        List<Double> gProST =  geneTreeProbability.calculateGTDistribution(currentST,geneTrees, null, false);
        return gProST;
    }

    //choose best gt for locus i
    public Double getBestGTi(List<Tree> gts, List<Tree> bestGTSi, Network tempST){
        Double ll = Double.NEGATIVE_INFINITY;
        List<Tree> tempGTS = new ArrayList<Tree>();
        for(int  i =0;i<gts.size();i++){
            Tree tempGT = Trees.readTree(gts.get(i).toNewick());
            tempGT.getRoot().setParentDistance(Double.NEGATIVE_INFINITY);
            tempGTS.add(tempGT);
        }
        //operator.scaleGTList(true, tempGTS,halfTheta);
        List<Double> llGTST = getGTSLLBySTYF(tempGTS,tempST);
        Tree best = gts.get(0);
        for(int i=0;i<gts.size();i++){
            Tree temp = gts.get(i);
            Double tLL = temp.getRoot().getParentDistance();
            if(ll<tLL*llGTST.get(i)){
                ll = tLL*llGTST.get(i);
            }
        }
        bestGTSi.add(best);
        return ll;
    }



    // P(Seq|gt) NOTICE: divergence time
    public double getSeqLLByGT(Tree gt, Alignment aln) {
        Utils._SUBSTITUTION_MODEL = "JC";
        Frequencies freq = new Frequencies(aln, false);
        //Trees.autoLabelNodes((MutableTree) gt.getTree());
        SubstitutionModel subst = new JukesCantor(freq);
        UltrametricTree ugt = new UltrametricTree(gt);
        BeagleTreeLikelihood beagle = new BeagleTreeLikelihood(aln, ugt, new SiteModel(subst), null);
        return beagle.calculateLogP();

    }

    //IMPORTANT: ST : tap -> cUnit -> tao    ||  return GTS: tao
    public List<Tree> simGTSByMS(Network<NetNodeInfo> currentST, double halfTheta, int numGTs) {
        Iterator<NetNode<NetNodeInfo>> it = currentST.bfs().iterator();
        List<Tree> gts = simulator.simulateGeneTrees(currentST, null, numGTs);
        return gts;
    }


}