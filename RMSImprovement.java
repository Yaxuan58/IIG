package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

/**
 * Created by doriswang on 6/26/17.
 * RAxML + ASTRAL
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
        _seqLens = new int[]{100};
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
                trueST = operator.load17Data("/Users/doriswang/PhyloNet/Data/17-taxon/004/ST0/target/", _seqLens[0], lociNum, trueSeq, trueGTS);
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
        System.out.println("\n" + "#0"  + " iteration start! ");
        System.out.println("Current LL is"  + GLOBALBESTLL);

        System.out.println("RF of Best ST : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)));
        System.out.println("Average distance of GTS is " + gtDis[0]);
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
        System.out.println("Best #iteration : " + bestLLIter );

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
        boolean changeST = false;
        if (STRATEGY == "NONE") {
            move = true;
        } else if (STRATEGY == "IMPROVE") {
            if (iter_LL >= maxLL) {
                move = true;
                changeST = true;
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
            if (iter_LL >= maxLL) {
                move = true;
                if (iter_LL > GLOBALBESTLL) {
                    maxLL = iter_LL;
                    changeST = true;
                    bestLLIter = itNum;
                    bestSTLL.add(iter_LL);
                    bestSTD.add(operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)));
                    bestITNum.add(itNum);
                    GLOBALBESTLL = iter_LL;
                    GLOBALBESTST = tempST.toString();
                    Iterator it = gts.iterator();
                    GLOBALBESTGTS = new ArrayList<String>();
                    while (it.hasNext()) {
                        GLOBALBESTGTS.add(it.next().toString());
                    }
                    System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GLOBALBEST ST will change");
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
//        if(itNum==1)
//            move = true;

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
//            if(changeST) {
//                GLOBALBESTST = tempST.toString();
//                System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GLOBALBEST ST changes");
//                System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ distance : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)));
//            }
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
        operator.getInitTree(trees,currentHGTLL,lociNum);
//        trees[0] =   "((14:0.31454687558361765776,((11:0.82162907848227140217,((2:0.21213617549968616371,(16:0.30874285527661915474,((O:2.16411663107923901705,5:0.00000100000050002909):0.13874749053706425528,7:0.12885641700209179983):0.04195215978636428200):0.03411192681770849194):0.09278690547910842246,(3:0.22773223241529208516,6:0.30382416641750692454):0.07121799895552061166):0.35993880396171901159):0.08436391279748052530,(10:0.29851327797909604778,((15:0.00699489578969816547,4:0.06878866163236331543):0.33242741866731395639,(8:0.04887987074296101103,(9:0.03934661035708268001,13:0.04565461459482478279):0.08552364758851636883):0.45532498856880831983):0.00000100000050002909):0.17969136349875336478):0.59944534893008827314):0.24999567242322306981,12:0.03760239523541118040,1:0.07470738375690583188):0.0;";
//
//        trees[1] =   "((12:0.12361676407710606163,(14:0.29854530790212269187,((11:0.70938783310738995347,(((16:0.20645685045080108821,(5:0.13111548320096755305,7:0.17786592620698446487):0.00000100000050002909):0.07596884052369075324,2:0.20586633980357621487):0.14488100386942084907,(6:0.22423614542996603283,3:0.33721160283012019399):0.09701882625037898389):0.31820826582494021739):0.19774848302319938953,(((15:0.05762723365359608690,4:0.00074397265823729136):0.21324453981056140650,10:0.24621417618781268888):0.20663788097302315805,(13:0.08853253815379567848,(8:0.24517296895700735337,9:0.11617089707114475194):0.00000100000050002909):0.40894511663032284687):0.41499072876581122316):0.23564586140236823408):0.28106915068299026528):0.01726442928107360664,O:2.51982782053959697421,1:0.00000100000050002909):0.0;";
//
//        trees[2] =   "(12:0.03278005088839869990,(14:0.56737459827483061048,(((3:0.22051508120226792631,6:0.32129173640134967016):0.10903118867408249926,((5:0.13996411818317891496,2:0.20949137493772729268):0.00680053087299913647,((16:0.07244104052657238113,((4:0.00000100000050002909,((10:0.40390219007271621621,(8:0.19492289773445739542,(13:0.06454070299038436320,9:0.05813845245107625292):0.08171098632263308625):0.64243101253526402328):0.07186234979492760189,15:0.03280472417950854036):0.01735731827216619935):0.99536400745019348868,O:1.17260912979677045520):0.48413008261810425159):0.20073138324380040975,7:0.17548509026967323843):0.05241231586506980844):0.07169381210951179806):0.42557212810089967148,11:0.78076711287401467931):0.21843936659289639457):0.43096273635074672237,1:0.00000100000050002909):0.0;";
//
//        trees[3] =    "(12:0.07814750913925462561,((((5:0.16393185638542814853,7:0.03540103824276032152):0.03875617974159424606,(((3:0.22740762980132450677,(O:3.12440103417684422737,6:0.00000100000050002909):0.24786085687006303990):0.17561568650961789650,(2:0.33118869684101687190,(10:0.04880558224466812617,((8:0.12811012465658619397,(9:0.09375474204761735653,13:0.03973135677494383372):0.08301802436893147841):0.47293492355626487456,(15:0.01705673241012171937,4:0.00741223453972161047):0.35290420064787214605):0.11394742041365805985):0.91955199824683342413):0.07618627936147999435):0.05010344146790666336,16:0.18675129130803563848):0.07540706174066456879):0.53068226649278515961,11:0.51109908338171028408):0.64659831110358112483,14:0.40657406969982035072):0.08188759754407447244,1:0.11590857187240341530):0.0;";
//
//        trees[4] =     "(((((O:1.38030592103746774235,((13:0.08870225792257502295,9:0.00518159941941468121):0.12675828544424624922,(8:0.17531367331599245762,(10:0.23463216458118990593,(4:0.04937977013478178634,15:0.00000100000050002909):0.20563657022660561391):0.19595329891032559066):0.09038860662259842749):0.09851505682910878525):0.37158746182814145209,((3:0.16588857651932034654,6:0.21607177341980249308):0.20386934608258866497,(16:0.24647979224018115585,(7:0.19865653832417351499,(2:0.29719239637448841007,5:0.21615008749966160839):0.02517531026384588200):0.06841706892048518562):0.02771452273417894518):0.38527135347237578822):0.26530444357339766359,11:0.59086736757480096127):0.54834626123268814801,14:0.44900367142044822488):0.15546876828543279569,12:0.06553336791484647150,1:0.03785019458015550370):0.0;";

//        currentHGTLL[0] = -2042.893488;
//        currentHGTLL[1] = -2090.553245;
//        currentHGTLL[2] = -2083.194595;
//        currentHGTLL[3] = -2057.147040;
//        currentHGTLL[4] = -2048.066367;

//        currentHGTLL[0] = -2084.536682;
//        currentHGTLL[1] = -2042.893488;
//        currentHGTLL[2] = -2090.553245;
//        currentHGTLL[3] = -2158.541632;
//        currentHGTLL[4] = -2035.803485;
//        currentHGTLL[5] = -2132.027419;
//        currentHGTLL[6] = -2083.194595;
//        currentHGTLL[7] = -2061.227966;
//        currentHGTLL[8] = -2057.147040;
//        currentHGTLL[9] = -2048.066367;

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
            System.out.println("# " + k + " : " + d);
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