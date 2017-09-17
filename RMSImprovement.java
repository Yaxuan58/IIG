package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

/**
 * Created by doriswang on 6/26/17.
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

public class RMSImprovement {
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
    private static int[] _seqLens;
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
    RMSImprovement() throws IOException, InterruptedException, ParseException {

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
        DATASET = "17";
        halfTheta = 0.04;//need to init in ALL classes ||  theta -> POP_SIZE
        lociNum = 10;
        ifOutGroup = false;

        ogHeight = new double[2]; // og height from all initGT[min, max]
        refineGTSize = 10;
        _scales = new double[]{1.0};
        _seqLens = new int[]{1000};
        taxaNum = 16;
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
        STRATEGY = "IMPROVE"; //RANDOM IMPROVE N
        //TODO Simulation
        simulator = new IIGTSimulator(lociNum, _scales, _seqLens, halfTheta, ITERATION);
        operator = new InferOperator(ITERATION);
        String resultFolder = RESULT_DIR + ITERATION + "/" + taxaNum + "_" + lociNum + "/" + refineGTSize + "/";
        operator.isExitsPath(resultFolder);
        if (ISREALDATA) {
            if (DATASET.equals(("17"))) {
                trueST = operator.load17Data("/Users/doriswang/PhyloNet/Data/17-taxon/004/ST3/1/", _seqLens[0], lociNum, trueSeq, trueGTS);
            }
        } else {
            List<String> trees = simulator.simulateData();
            trueST = trees.get(0);
            trueGTS = trees.subList(1, trees.size());
            for (int i = 0; i < lociNum; i++) {
                Alignment aln = operator.loadLocus(i, _seqLens[0], taxaNum);
                trueSeq.add(aln);
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
        Tree inferredST = operator.inferSTByGLASS(gts);
        if(inferredST.getNode("1").getParentDistance()!=Double.NEGATIVE_INFINITY)
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
//            if(inferredST.getNode("1").getParentDistance()!=Double.NEGATIVE_INFINITY)
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
        operator.getUpdateSh(lociNum,refineGTSize);
        String trees[]  = new String[lociNum];


//
//        trees[0] = "(15:0.13667690038143084919,(6:0.17941588735625457751,(((13:0.18314307288981809818,(14:0.10470543838806953274,(O:2.52983082977044659856,9:0.00000100000050002909):0.09919410343836373989):0.02667551081556799111):0.59822222844529737706,(7:0.32001398458337937392,(5:0.19700960483974505610,2:0.13889895039793248577):0.24072281015341878696):0.54045228728200411794):0.47501940067007164537,((((12:0.06457479862596982034,8:0.05047967511516762706):0.01528964148654596983,11:0.08857283882026469046):0.12285855308910910433,(16:0.04579972051457772281,4:0.04977829012286350691):0.13403294260698458973):0.04916680670533749020,(3:0.03699036101362166568,10:0.07025923800845568223):0.19571900808962378049):0.04935131505839653737):0.15384066673559251837):0.09137427125469683997,1:0.11850332944440154781):0.0;";
//
//        trees[1] =   "((((3:0.07456131619974389058,10:0.08686088469544360480):0.11543050060333480411,6:0.21337555813865852206):0.03593805565211598968,(((12:0.07619873619562866907,8:0.07330221679152214287):0.02423034390024219123,((16:0.05448361662574916636,4:0.04500088030029123637):0.02971969101069311253,11:0.08427785218204174311):0.02350929415679824455):0.26400839651441326827,(((5:0.19353910530420978531,2:0.17267677433425410172):0.25393165391092359373,(O:2.59490630106479125416,7:0.00000100000050002909):0.42871908254518936721):0.50568942316849074814,(14:0.06487550697360448815,(9:0.10390279644522203151,13:0.24959979822419040008):0.02604750387481960180):0.65901794972666727368):0.57485056656897948191):0.03432168180356610226):0.21840656437980787707,15:0.01758637479059742964,1:0.01410506451840563488):0.0;";
//
//        trees[2] =   "((((((13:0.51976386343081337760,(14:0.13330180945183303654,9:0.12787592100154090979):0.34489186529618587329):0.24389962845940790781,((2:0.17629597087119192045,5:0.17973420467297737679):0.20336164025509062547,7:0.35886959759736403175):0.54932050023449863208):0.60186478401648701997,(((16:0.08859877448593381222,(11:0.08576391602261243530,8:0.08725356983964804858):0.01323589645132520511):0.00292725175599889333,4:0.09480666753383724943):0.01585264426388520442,12:0.12065365382775451730):0.18183348364197154945):0.04143020041317303781,((O:2.88791842696944067725,10:0.00000100000050002909):0.15264838547301565197,3:0.09147219223273143907):0.09256240342943281019):0.01531658694594161332,6:0.21459220588999797230):0.15510582587866750259,15:0.05121892370926374449,1:0.05169090246069189820):0.0;";
//
//        trees[3] =    "((((3:0.03328478508260293567,10:0.03390907345615705448):0.20259843116926987139,((((13:0.12491425830734259605,(9:0.13224310587806728523,14:0.15880350018408262436):0.06516668220557707281):0.66114858482079852386,(7:0.46545180507966676942,(5:0.17286594599072691025,2:0.20798340985116228752):0.25171098799788349520):0.36092326313796624371):0.81974063132489194405,(12:0.12947822452009483118,(((16:0.07455251235318732794,4:0.06167040738485125156):0.01775944829558460253,11:0.07730165063638035072):0.00729118140273579338,8:0.09499998819827583374):0.03667998503775647678):0.16235142500039784208):0.06225845451056837704,6:0.23552132008960657839):0.00173114907783549561):0.16868939226373022566,15:0.03446359353165107803):0.03381872052349579882,O:2.59110643698017284464,1:0.00000100000050002909):0.0;";
//
//        trees[4] =     "(15:0.01943383450403767670,((((11:0.06179882803883740561,(12:0.05035750198205708739,8:0.07509567994214405573):0.00922125810955304499):0.01614637813399492916,(16:0.04898703377623826738,4:0.04032381466130432329):0.03584105916704791867):0.18770870341224918376,((7:0.44342015007137436378,(5:0.22274280142839283925,2:0.23252453466323363962):0.23287132388268302896):0.48961791314592278956,((O:2.83546277797208201576,14:0.02060267341399979121):0.22034958539052562632,(9:0.16348532865513493162,13:0.14682900035577070574):0.04036225719565540715):0.38364064026788158301):0.54981051763459654058):0.03040075619640988241,((10:0.10846557456973213163,3:0.12322296304605627471):0.08093511899902580686,6:0.15594752744731857352):0.03822583811945361870):0.21801955664206656982,1:0.01543210187050158416):0.0;";
//
//        trees[5] = "(15:0.00397082009400540666,((((2:0.25448483315771308622,(7:0.37304784410528429861,5:0.27039221170036381592):0.01843634299435019039):0.37731660391483623496,(13:0.13797755170391745594,(14:0.15772310599215921822,9:0.14151154893142200986):0.08531220735721462001):0.64887084594286725814):0.53048598889066256312,(11:0.20259172875819372872,((12:0.04522955153799452371,(O:3.62437759798986336435,8:0.00000100000050002909):0.03456392459083565405):0.06543968908926288408,(16:0.06132262133093489831,4:0.05899374354053352199):0.04904779038469816305):0.05171146455838783462):0.15041073292211745716):0.10471251052268386317,((10:0.14867873462266154028,3:0.19179131856235942521):0.02062958586818367315,6:0.21830122942385349338):0.00012462281084881831):0.20497471458920793475,1:0.00003800775703741025):0.0;";
//
//        trees[6] =   "(15:0.01340719788672756034,((3:0.11702866686351730274,10:0.12105421267297565568):0.11274277958943153266,(6:0.23687079260408075432,(((14:0.11778549848637165365,(9:0.14033462745518815451,(O:2.76249707226012475303,13:0.00000100000050002909):0.12775347422016711252):0.05619923850222211664):0.60215345351389981676,(7:0.62953911084534286413,(2:0.20997641999044486250,5:0.19762035305007619179):0.28166564975198654652):0.43158640202017617415):0.55991397193039593727,(11:0.08508033326340502800,((4:0.05274640748364878978,16:0.05453802296561190144):0.05125446592287823122,(12:0.07963758317124204622,8:0.08948006398259898686):0.00605957326175058774):0.02060750863767717156):0.16053370787286891375):0.09000065158386290343):0.01909507939924651176):0.21153902364550911108,1:0.00988132294443915676):0.0;";
//
//        trees[7] =   "(15:0.00148942483319117192,(((3:0.07901037708243154045,10:0.07538999499584302100):0.16909044882311294322,6:0.24704197431091723391):0.02203172295185615998,(((11:0.06733559329529512083,8:0.07034527231192243268):0.02457055391727401561,(12:0.08490935533002856583,(4:0.03488140446756029384,16:0.02751370936714118628):0.05855662724819503395):0.00035729501195154158):0.14748928992496385115,(((14:0.12950891828806404549,9:0.09420933365786403335):0.26342607481691898696,13:0.32422974716831293174):0.42514923235006252256,((O:2.64783515442788264238,7:0.00000100000050002909):0.25489699233573032133,(2:0.17728324222978675273,5:0.18974887366867471128):0.19337733522390929131):0.61875156253788388128):0.48101668914271300537):0.02446458206623858267):0.24815267201201385294,1:0.00252359363361547495):0.0;";
//
//        trees[8] =    "((((((8:0.07284344508783495431,11:0.06747386570266492023):0.01153391189407957165,12:0.08731406446495487694):0.08201328106007992058,(16:0.07609769574691246929,4:0.06326240185955851747):0.07882629918881599307):0.12632844888304320397,((9:0.12593890620588843454,(14:0.15972521544445245079,13:0.16568705756847376676):0.02704581125835135713):0.55112868044337770534,((5:0.32399373243079776508,(O:34.53877639491068407551,2:0.13954058482321829437):0.13536553160197337120):0.03976625268769121840,7:0.38648787856554217557):0.60147607181151452327):0.52793707539480827506):0.10778715443977426869,(6:0.22312684207767938260,(10:0.04155914040006220000,3:0.02995864690887481072):0.17518705093982492316):0.05102018892211454548):0.21307375585484755742,15:0.00679916994438379153,1:0.00630747249608694643):0.0;";
//
//        trees[9] =     "(15:0.03333102042810749677,(6:0.21207556760119392059,((((12:0.09436627552704630884,(8:0.09760191203130563764,4:0.08475227802637472219):0.00356048590444453441):0.01604237590230496099,(16:0.10009239971457528640,11:0.11702952280324825651):0.00846070571061136034):0.22971039519346816982,(((5:0.19435686247227551560,2:0.16777469240444717324):0.20903491081405947694,(O:2.83158114621593748339,7:0.00000100000050002909):0.41096991831639123616):0.46638545808145021043,((9:0.12640993484282339598,14:0.12318870773738074820):0.03561002125656068107,13:0.15282685905869428100):0.62207246045142272628):0.54943381491707543596):0.04221726061833907390,(10:0.19598269074554072744,3:0.20310452098175915681):0.05238864995639432726):0.02886720575213044468):0.16708794324032599787,1:0.02510276907824401244):0.0;";



        trees[0] = "(((8:0.12665072233407143054,(13:0.07874862553094696960,(4:0.06553391471435220106,16:0.06118176905784743908):0.00336624455713347521):0.05561757777289377813):0.08370933109144682027,((((11:0.12538772589814314862,12:0.14729176352543021933):0.19925818747059068081,(6:0.28266395897161594730,(14:0.26784969431235500226,(9:0.15678911972206122094,((7:0.04914790266946252040,15:0.03857135956628036932):0.05714315701900665734,3:0.11465684206308847382):0.09294054220211653283):0.06515039034222391834):0.07734907311862002088):0.31717514959328635404):0.69781466906033051689,(10:0.00000100000050002909,O:2.54875493548852460890):0.51568739479809344584):0.34115909979772196525,5:0.42532518977683431149):0.40661716912768269827):0.10965767299309953087,2:0.03029222628298596995,1:0.03090745685304375351):0.0;";

        trees[1] =   "(((4:0.07161737792867628460,16:0.07814557358898857609):0.01197053724567953385,(O:2.97820055153061113984,13:0.00000100000050002909):0.07251567395339961197):0.04682244136314300115,(((((11:0.17045874179272871296,12:0.16575489937724716127):0.30668149964041974131,(6:0.35807070825338144227,(14:0.23557099149862559262,(9:0.20711550153300420374,(3:0.16588375150455242357,(7:0.04397216789025047701,15:0.05580425427373424119):0.11092205382868802754):0.02856045210653231775):0.07878684044102474138):0.07589550989989983509):0.06539725373422251853):0.68290462788643779835,5:0.50092433967929228533):0.09269360180416084949,10:0.43502123451005125787):0.22244533870530724950,(8:0.10180104805249830247,2:0.14198964324861595987):0.06032964739967637546):0.02560351188369941466,1:0.13150902214894305553):0.0;";

        trees[2] =   "(2:0.01268885153345708332,(((16:0.06536451733972177547,13:0.05617218203853491937):0.06957633331599036297,(((((11:0.23754759528454022499,12:0.25617740747540557189):0.21483359523667738822,(6:0.27908574211599967985,(14:0.24987601310472482008,(9:0.17774260898159230382,(3:0.10964937066662706899,(7:0.05391299799707125612,15:0.05631952208833961721):0.07350687439300532255):0.08371962079551839575):0.03987934426846617697):0.07473528769521363413):0.17947532498575516358):0.66747164429602767743,5:0.46557147426482464736):0.16938674446385451611,10:0.41765387789289787346):0.26540867588744704308,(O:2.61617871811515190217,8:0.00000100000050002909):0.13758142951664192877):0.00720466809301000856):0.03222462417951449387,4:0.09863931840416373387):0.10543356449880253445,1:0.01168803035880548562):0.0;";

        trees[3] =    "(2:0.02297155276733786611,((((5:0.44900814112594805705,((((((15:0.02626772505762545520,7:0.03189989830739448778):0.08614512304840943324,3:0.13369380116305951134):0.07962089325932954209,9:0.18711424243231641884):0.09174046481531067387,6:0.32323058176349528381):0.04808537743644625273,14:0.39404510057353597796):0.16564668572872970342,(11:0.16132406644932253958,12:0.15543635298376717868):0.40360710602251143042):0.80094467798724966467):0.20207485912397649130,10:0.45423923187838560755):0.27923713177475789626,(16:0.04622041194169525247,(4:0.05590218706146036359,(O:34.53877639491068407551,13:0.03073667873133120038):0.03063629087975586196):0.01336091164991447082):0.05257587883829752118):0.02823861560184830893,8:0.09119664558494343376):0.08247624823953297879,1:0.02460938139618533069):0.0;";

        trees[4] =     "(2:0.00895547174244147592,((((5:0.57599515365863418470,((((((15:0.02784349450310677140,7:0.02073491155114557186):0.12304833838322334905,3:0.12817321943102549797):0.10656525562517121797,9:0.21566785928155485230):0.04459373119931073909,14:0.32442061063846694413):0.03929039004877590835,6:0.27010012550325995795):0.13624024662695782828,(11:0.16298369047824012856,12:0.15078149662412818066):0.24177260817645124247):0.76313872155356077265):0.17584820698157610508,(O:2.38292802939331016532,10:0.00000100000050002909):0.42680035700721835612):0.09171430281495489278,(4:0.07712772992455636800,(16:0.06486056922625044185,13:0.07706349754631856841):0.01006282801576144185):0.17675776328517298075):0.13598742712415126666,8:0.14511879047540379717):0.14626796756720114590,1:0.02379731859851487441):0.0;";

        trees[5] = "(2:0.00716741653811027793,(((O:2.81063139513574800077,13:0.00000100000050002909):0.07153445345682026002,(16:0.07579278770047610803,4:0.05041810912628691410):0.00426562427380143376):0.10321954058821627975,(((5:0.48744261234514218994,((11:0.16389690428297129521,12:0.15930972017128952611):0.31241853418647347329,((14:0.25412066644117686787,(9:0.18118296528665517320,(3:0.16195962385367401737,(15:0.02510481566882414792,7:0.01502168978107644610):0.13565415839932379316):0.04346623023238371797):0.07142322276574133988):0.13803180012244167796,6:0.42766052210111948906):0.11925155208606230683):0.61321597593356202971):0.15904829696769368086,10:0.35443263750672981960):0.20372155075369946520,8:0.20013603641612159878):0.02918443144997560912):0.11966842656841403247,1:0.00389773719092853304):0.0;";

        trees[6] =   "(2:0.06869092212203486980,(((16:0.11739061147114622508,(4:0.09691420584394899729,13:0.09788062279101478191):0.00563016663913393049):0.02473537602687249717,((5:0.34967039052202569804,((((((15:0.08693671734387813099,7:0.08002216896969591153):0.11144317457913843761,3:0.17197925175882633475):0.02636586654117735115,9:0.18159682796863721443):0.08758180140866186902,(O:2.88368583952718537944,14:0.00000100000050002909):0.34134547627472788811):0.03497564908045690257,6:0.33739774111990006311):0.18367575470966121531,(11:0.18133746297583031026,12:0.15329631611469227148):0.24125089798325854984):0.68855338295960810679):0.18703063050920445964,10:0.40073105863775759783):0.28967593024626420162):0.03148436888006630324,8:0.08232338869505423373):0.03343366450991742261,1:0.04946479871327604666):0.0;";

        trees[7] =   "(2:0.05252554541672248972,(8:0.10257225852806002575,((4:0.07963287986488684633,(16:0.03351635914294477658,(O:2.66398103264252394951,13:0.00000100000050002909):0.02382872048571542042):0.07208533080845626329):0.01724809124727976969,((((11:0.22633061735356621580,12:0.18910677142294668451):0.26744804396034854710,(6:0.29886517321447148010,(14:0.22963151559368591426,(9:0.18408241571828720673,(3:0.12325136258092965480,(7:0.07454925247804429500,15:0.08288643490664007740):0.03165017609223371581):0.07806207283952658338):0.05160169076451684433):0.06459881618583092344):0.12893789962714469199):0.63342786123994332392,5:0.52919883118548649570):0.16399995535176439820,10:0.44030021380866857239):0.23952878155451978737):0.00615439314224362072):0.08744822077699568319,1:0.04633357207094759844):0.0;";

        trees[8] =    "((8:0.09357975105865952437,((4:0.06986192983561885506,(13:0.01415678481832776010,16:0.02778896180588393341):0.03449734335342628105):0.02279061618853505181,((5:0.54855913524371957557,((11:0.17024134152140074638,12:0.18754816881268715489):0.26650551909082798030,(6:0.37361183473615677819,((O:4.15248019894930830986,14:0.00000100000050002909):0.28986912521193031189,(9:0.20632043504415326085,(7:0.12228579100754796705,(15:0.08840061034428603470,3:0.10376984570858960311):0.03374281090549418466):0.07467203630689905891):0.09829670685448234357):0.04670227790615809338):0.05472892551987983134):0.53413029202221029568):0.14643376746244468611,10:0.46733281444325486254):0.32212622855880734818):0.00914347498186208509):0.07270220052379439135,2:0.02660116675591434401,1:0.01445213941818990037):0.0;";

        trees[9] =     "((2:0.02814690778669355284,(((4:0.06136954254347021243,(16:0.03206619047682021661,13:0.04300620435389011786):0.02385623149958504899):0.17927366497153376934,((((11:0.17073958075256395750,12:0.18428095547334832904):0.29563112590624296416,(6:0.30585680010332488843,(14:0.29618490757751658737,(9:0.20069729616466130961,(3:0.10425738266804403243,(7:0.02156414928427469874,15:0.03064841900677113193):0.07433896413031769457):0.08149193985525619854):0.03746734169704837858):0.12945738577505105926):0.04916906106607850380):0.86438955069204337356,5:0.71350430374192108651):0.08491518464249528386,10:0.41047454637887142992):0.21002297172168227224):0.09729492694935634733,8:0.11744915279944297126):0.08200942657062097740):0.02340815097691369140,O:2.69990631054782337372,1:0.00000100000050002909):0.0;";
        List tempGTs = new ArrayList<Tree>();
        double ogMin = Double.POSITIVE_INFINITY;
        double ogMax = 0.0;
        for(int i=0;i<lociNum;i++){
            Tree tempT = Trees.readTree(trees[i]);
            tempT = operator.rerootRAxML(tempT,"O");
            //Tree tempT = tempOGTs.get(i);
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
        //operator.scaleGTList(false, tempGTs, halfTheta);
//        currentHGTLL[0] = -14957.630015;
//        currentHGTLL[1] = -15238.757075;
//        currentHGTLL[2] = -15143.478340;
//        currentHGTLL[3] = -14499.177003;
//        currentHGTLL[4] = -14602.864392;
//        currentHGTLL[5] = -14669.564935;
//        currentHGTLL[6] = -14539.299872;
//        currentHGTLL[7] = -14887.092847;
//        currentHGTLL[8] = -14644.965958;
//        currentHGTLL[9] = -15246.248191;

        currentHGTLL[0] = -15050.871093;
        currentHGTLL[1] = -15489.769359;
        currentHGTLL[2] = -15021.499400;
        currentHGTLL[3] = -14758.680813;
        currentHGTLL[4] = -15316.116644;
        currentHGTLL[5] = -14957.244398;
        currentHGTLL[6] = -15413.22865;
        currentHGTLL[7] = -14860.775460;
        currentHGTLL[8] = -14602.886459;
        currentHGTLL[9] = -14908.402736;
//        currentHGTLL[0] = -17953.316158;
//        currentHGTLL[1] = -18116.952650;
//        currentHGTLL[2] = -17866.466959;
//        currentHGTLL[3] = -17372.748495;
//        currentHGTLL[4] = -17811.657316;
//        currentHGTLL[5] = -17718.722821;
//        currentHGTLL[6] = -17854.758719;
//        currentHGTLL[7] = -17531.358452;
//        currentHGTLL[8] = -17396.247747;
//        currentHGTLL[9] = -17737.643832;
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
        trueST = operator.rerootAndRemove(trueST,"O");
        stDis[0] = operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
        double thisD = 0.0;
        //finalGTS = tempGTs;
        for (int k = 0; k < lociNum; k++){
            trueGTS.set(k,operator.rerootAndRemove(trueGTS.get(k),"O"));
            double d = operator.getDistance((Tree) finalGTS.get(k), Trees.readTree(trueGTS.get(k)));
            gtDis_Locus[0][k] = d;
            msDist_Locus[0][k] = 0;
            thisD += d;
        }
        gtDis[0] = thisD / lociNum;
        return INIT_ST;
    }


    public static void main(String[] args) throws IOException, ParseException, InterruptedException {

        RMSImprovement rms = new RMSImprovement();
        InferOperator operator = new InferOperator(ITERATION);
        IIGTSimulator simulator1 = new IIGTSimulator(2, _scales, _seqLens, 0.04, 5);
        //ast.initST(trueSeq);
        rms.iigt(trueSeq,ITERATION);
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