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

public class ASTImprovement {
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
    private static List<Alignment> true16Seq;

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
    private int stIndex; //index of true gt from /Users/doriswang/PhyloNet/Data/17-taxon/

    //All parameters should be initialed here
    ASTImprovement(int index) throws IOException, InterruptedException, ParseException {

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
        ITERATION = 200; //{1k, 5k ,10k}
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
        stIndex = index;
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
            if (DATASET.equals(("17"))) {
                trueST = operator.load17Data("/Users/doriswang/PhyloNet/Data/17-taxon/004/ST0/1/", _seqLens[0], lociNum, trueSeq, trueGTS);
                true16Seq = trueSeq;
                //true16Seq.remove()
            }
            if (DATASET.equals(("RANDOM"))) {
                //trueST = operator.loadRandomData(stIndex, _seqLens[0], lociNum, trueSeq, trueGTS);
            }
        } else {
            List<String> trees = simulator.simulateData();
            trueST = trees.get(0);
            trueGTS = trees.subList(1, trees.size());
            for (int i = 0; i < lociNum; i++) {
                Alignment aln = operator.loadLocus(i, _seqLens[0], taxaNum, "/Users/doriswang/PhyloNet/Data/17-taxon/004/ST0/1/Seq/");
                trueSeq.add(aln);
            }
        }
    }

    public static void main(String[] args) throws IOException, ParseException, InterruptedException {

        ASTImprovement ast = new ASTImprovement(0);
        InferOperator operator = new InferOperator(ITERATION);
        IIGTSimulator simulator1 = new IIGTSimulator(10, _scales, _seqLens, 0.04, 100);
        //ast.initST(trueSeq);
        ast.iigt(trueSeq,ITERATION);
        System.out.println("Finish");
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

        //cGTS: only use topo
        List<Tree> currentGTS = simulator.generateGTS(tempST,refineGTSize, _msdir,null, cHeights, false); //operator.MS
        operator.getMSDist(trueTopos, currentGTS, msDist_Locus, itNum);
        //List<String> topos = operator.getMSTopos(currentGTS, checkedTopo);

        List<Tree> gts = getBestGTSByMS(currentGTS,tempST);
        boolean move = false;
        double p2 = currentLL - currentHLL;
        currentITLL[itNum] = iter_LL;

        if (FULL_LL)
            iter_LL = currentHLL*ll_Ratio + p2*(2-ll_Ratio);
        else
            iter_LL = currentHLL;

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
                }
            } else {
                double random = Math.random();
                if (random < 0.975)
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
            Tree inferredST = operator.inferSTByGLASS(gts);
            if(inferredST.getNode("1").getParentDistance()!=Double.NEGATIVE_INFINITY)
                gts.add(inferredST);
            List<Tree> gtsOG = operator.getOGGTS(gts, ogHeight[0]);
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

        //List<Tree> tempOGTs = operator.initByRAxML(aln, _seqLens[0],taxaNum, lociNum);
        String trees[]  = new String[lociNum];
        trees[0] = "(15:0.06385239050018931550,((((7:0.35507858635141326120,(5:0.29353650672204489869,2:0.23406687701105596822):0.06643191841041275192):0.43443454037783985067,((14:0.12499816664032974145,9:0.14170502523795053262):0.04374612299215879102,13:0.15473405018419678081):0.52957691547077978544):0.53786789459448192119,(8:0.12080737028319389614,((12:0.08356589617662767144,11:0.06447326508961689906):0.03147560317215407355,((O:2.34057499021339276979,4:0.00000100000050002909):0.05174336931150425728,16:0.05053598520394290278):0.04770292865750219713):0.03704580009898181259):0.14173872698785303093):0.05301445258794363485,((3:0.06392832764165792825,10:0.05796107858634157867):0.12011281425939168699,6:0.20601493580392435390):0.07731084362708398605):0.24172597823533428785,1:0.08367374764129421882):0.0;";

        trees[1] =   "(((3:0.07067676949193306413,10:0.05490362740030942101):0.15915501529849843121,(6:0.18583036929393462189,((((12:0.03365158635464698955,8:0.05339983396435030866):0.07290109521790824609,11:0.08699423187957987247):0.30145138455352843421,((((O:2.79339209973832547362,5:0.00000100000050002909):0.32374386394702830971,2:0.30392381730752554558):0.04923391482203538783,7:0.38808719386447132482):0.41867808904179854013,((13:0.13693194675898351376,14:0.13127426326912855492):0.05582651902402614025,9:0.20796222415194326461):0.64009955891090419833):0.52676112470034985868):0.10992175270616143246,(16:0.07230998576624945995,4:0.06586836885251330653):0.23735830000761412584):0.05517271419130038296):0.02295560794657667755):0.23447979019970513481,15:0.00599934813918216042,1:0.00406576464609070994):0.0;";

        trees[2] =   "(((6:0.21864692152097200961,(((16:0.09948079467352581584,(12:0.06811430488209768708,8:0.07431607512919138903):0.02648238011667029476):0.03778550672047790993,(4:0.13060908038257934560,11:0.08562746962628120517):0.01219604300873776588):0.14581266225694869410,((13:0.16237848382081532250,(O:2.34077418772927670076,(14:0.14089121826408579352,9:0.11920005236362521406):0.16860184982863507530):0.02427074271366035096):0.59314355288001452404,((5:0.21582170498353936416,2:0.23004171136387729923):0.10360571271432542084,7:0.35052215593934077376):0.61248499378284471106):0.65216579736076796259):0.05585956028786604605):0.00668098040161886706,(3:0.14093840019188050294,10:0.14105455128389468578):0.06270524915762901164):0.15780604191116268109,15:0.04601464874506398101,1:0.03925786701673324858):0.0;";

        trees[3] =    "(15:0.01504073197183840652,((((((2:0.27715353055650876479,5:0.23994623184437588459):0.16960311258998658768,7:0.27429391383418738615):0.57488257589166602557,(13:0.09258402465310638929,(14:0.09665510942607179190,(O:3.03517750932576335643,9:0.00000100000050002909):0.11404251417786588629):0.09216691812612336221):0.56963286005863877914):0.43066569270239335454,(11:0.13523535475148720786,(12:0.09422117869190024808,((16:0.06566091067862804553,4:0.06232752953366831050):0.02686209599698902856,8:0.08478586141789691555):0.01379996860303076521):0.03761460486523916097):0.17379066189132019726):0.10499904924467018397,(3:0.04023156769678863653,10:0.03203571117735003887):0.23644419357652782243):0.01263779425747846991,6:0.22461533901501193333):0.21618297005840378389,1:0.01669062174999247361):0.0;";

        trees[4] =     "((((((9:0.14346674889858623825,14:0.15915614783964654455):0.05233336282342995749,13:0.17119692549367435119):0.45529474007969017446,((2:0.15830736928244024120,5:0.19584476185184726549):0.25736835478331748783,(O:2.48142909549519119494,7:0.00000100000050002909):0.36180239062652264082):0.61580560652742155625):0.52341327161783202104,(11:0.09770273434829912507,(12:0.09415447162603486631,(8:0.08649011051468645517,(4:0.06468901630280333992,16:0.07160461603597464975):0.02931889857258338991):0.00345164581426198741):0.01330590757513378138):0.22187787603728831298):0.12127200951536722118,((10:0.04019493585650776163,3:0.04866162700196443452):0.17646419724619608060,6:0.21046112688285797954):0.01347532607852260957):0.21839907362733285145,15:0.02087511746025366710,1:0.02122372485415865984):0.0;";

        trees[5] = "(15:0.06385239050018931550,((((7:0.35507858635141326120,(5:0.29353650672204489869,2:0.23406687701105596822):0.06643191841041275192):0.43443454037783985067,((14:0.12499816664032974145,9:0.14170502523795053262):0.04374612299215879102,13:0.15473405018419678081):0.52957691547077978544):0.53786789459448192119,(8:0.12080737028319389614,((12:0.08356589617662767144,11:0.06447326508961689906):0.03147560317215407355,((O:2.34057499021339276979,4:0.00000100000050002909):0.05174336931150425728,16:0.05053598520394290278):0.04770292865750219713):0.03704580009898181259):0.14173872698785303093):0.05301445258794363485,((3:0.06392832764165792825,10:0.05796107858634157867):0.12011281425939168699,6:0.20601493580392435390):0.07731084362708398605):0.24172597823533428785,1:0.08367374764129421882):0.0;";

        trees[6] =   "(((3:0.07067676949193306413,10:0.05490362740030942101):0.15915501529849843121,(6:0.18583036929393462189,((((12:0.03365158635464698955,8:0.05339983396435030866):0.07290109521790824609,11:0.08699423187957987247):0.30145138455352843421,((((O:2.79339209973832547362,5:0.00000100000050002909):0.32374386394702830971,2:0.30392381730752554558):0.04923391482203538783,7:0.38808719386447132482):0.41867808904179854013,((13:0.13693194675898351376,14:0.13127426326912855492):0.05582651902402614025,9:0.20796222415194326461):0.64009955891090419833):0.52676112470034985868):0.10992175270616143246,(16:0.07230998576624945995,4:0.06586836885251330653):0.23735830000761412584):0.05517271419130038296):0.02295560794657667755):0.23447979019970513481,15:0.00599934813918216042,1:0.00406576464609070994):0.0;";

        trees[7] =   "(((6:0.21864692152097200961,(((16:0.09948079467352581584,(12:0.06811430488209768708,8:0.07431607512919138903):0.02648238011667029476):0.03778550672047790993,(4:0.13060908038257934560,11:0.08562746962628120517):0.01219604300873776588):0.14581266225694869410,((13:0.16237848382081532250,(O:2.34077418772927670076,(14:0.14089121826408579352,9:0.11920005236362521406):0.16860184982863507530):0.02427074271366035096):0.59314355288001452404,((5:0.21582170498353936416,2:0.23004171136387729923):0.10360571271432542084,7:0.35052215593934077376):0.61248499378284471106):0.65216579736076796259):0.05585956028786604605):0.00668098040161886706,(3:0.14093840019188050294,10:0.14105455128389468578):0.06270524915762901164):0.15780604191116268109,15:0.04601464874506398101,1:0.03925786701673324858):0.0;";

        trees[8] =    "(15:0.01504073197183840652,((((((2:0.27715353055650876479,5:0.23994623184437588459):0.16960311258998658768,7:0.27429391383418738615):0.57488257589166602557,(13:0.09258402465310638929,(14:0.09665510942607179190,(O:3.03517750932576335643,9:0.00000100000050002909):0.11404251417786588629):0.09216691812612336221):0.56963286005863877914):0.43066569270239335454,(11:0.13523535475148720786,(12:0.09422117869190024808,((16:0.06566091067862804553,4:0.06232752953366831050):0.02686209599698902856,8:0.08478586141789691555):0.01379996860303076521):0.03761460486523916097):0.17379066189132019726):0.10499904924467018397,(3:0.04023156769678863653,10:0.03203571117735003887):0.23644419357652782243):0.01263779425747846991,6:0.22461533901501193333):0.21618297005840378389,1:0.01669062174999247361):0.0;";

        trees[9] =     "((((((9:0.14346674889858623825,14:0.15915614783964654455):0.05233336282342995749,13:0.17119692549367435119):0.45529474007969017446,((2:0.15830736928244024120,5:0.19584476185184726549):0.25736835478331748783,(O:2.48142909549519119494,7:0.00000100000050002909):0.36180239062652264082):0.61580560652742155625):0.52341327161783202104,(11:0.09770273434829912507,(12:0.09415447162603486631,(8:0.08649011051468645517,(4:0.06468901630280333992,16:0.07160461603597464975):0.02931889857258338991):0.00345164581426198741):0.01330590757513378138):0.22187787603728831298):0.12127200951536722118,((10:0.04019493585650776163,3:0.04866162700196443452):0.17646419724619608060,6:0.21046112688285797954):0.01347532607852260957):0.21839907362733285145,15:0.02087511746025366710,1:0.02122372485415865984):0.0;";

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

        //operator.scaleGTList(false, tempGTs, halfTheta);
        currentHGTLL[0] = -14957.630015;
        currentHGTLL[1] = -15238.757075;
        currentHGTLL[2] = -15143.478340;
        currentHGTLL[3] = -14499.177003;
        currentHGTLL[4] = -14602.864392;
        currentHGTLL[5] = -14669.564935;
        currentHGTLL[6] = -14539.299872;
        currentHGTLL[7] = -14887.092847;
        currentHGTLL[8] = -14644.965958;
        currentHGTLL[9] = -15246.248191;
//        currentHGTLL[0] = -14885.327181;
//        currentHGTLL[1] = -14376.336898;
//        currentHGTLL[2] = -15231.805277;
//        currentHGTLL[3] = -14575.994606;
//        currentHGTLL[4] = -14400.796516;
//        currentHGTLL[5] = -14776.970112;
//        currentHGTLL[6] = -14702.383736;
//        currentHGTLL[7] = -14283.969575;
//        currentHGTLL[8] = -14642.651198;
//        currentHGTLL[9] = -14987.974635;


        //        trees[0] = "((((11:0.70644398931864882396,(3:0.34767924403943067935,(((5:0.15710353984118755055,7:0.18935903425796674071):0.08874855967111235089,(16:0.20087264366883630928,2:0.25783647194418191084):0.03843178499567419248):0.05896998829160856465,6:0.38085907307506022823):0.07311579398713173583):0.37200568204787165527):0.01921612205860232306,(((13:0.07504554726485214433,9:0.06313470059355173747):0.06576799914386330381,8:0.14254521883365151025):0.26792964654666351354,(10:0.25828580936016115599,(15:0.01424732461309747601,4:0.01322523820803084803):0.18785549914016130835):0.19246845125134842691):0.55215562996529365947):0.46131833154437729139,(O:4.05911231085117574224,14:0.00000100000050002909):0.29596406113180256936):0.28137585062420888971,12:0.01731115410631995505,1:0.02501476357694640967):0.0;";
//
//        trees[1] =   "(12:0.04399827182659500624,(((((10:0.25444716531487121536,(15:0.00355401690007028448,4:0.00549988354214483999):0.19366689258000677198):0.18235892591196123380,(8:0.11046727708066554630,(9:0.07488860481588313633,13:0.04400698866155940830):0.07220305126224457748):0.22117053026417671235):0.54084844484732019776,((2:0.25676556371336750617,(16:0.19503163588057048639,(5:0.14646252764275405345,7:0.19104839503609721829):0.08411985892033234380):0.06458311226381635572):0.07172982922310700005,(3:0.30548101010256234478,6:0.31611224179943486767):0.03510270822114055361):0.36270266518108468246):0.00000100000050002909,(O:2.77115580175478282143,11:0.00000100000050002909):0.62484388145985103602):0.43031475524696627399,14:0.25858177305660351442):0.27488441331466606510,1:0.03821119722230179039):0.0;";
//
//        trees[2] =   "(12:0.05389264015061668439,(14:0.31036431233883499115,((((8:0.12617679588613076036,(13:0.05779635336035422732,(O:2.36442923010710437381,9:0.00000100000050002909):0.05984303062817749458):0.07391840543060251845):0.34912355218082097297,(10:0.20916395904115100968,(4:0.02087354149109655471,15:0.03389591261305226894):0.16063257781127360246):0.19088129202874004609):0.48533178418738748272,((2:0.19709442058381326790,(16:0.22116446851592569400,(7:0.14106565037434873333,5:0.14933792283054186933):0.09067034675596909177):0.02074605264692864354):0.20315483363510919435,(6:0.17929711217024171810,3:0.22547757042749946854):0.03814184753746297646):0.54190226689483522726):0.06556032698135033643,11:0.74641835382036614277):0.34983134815328109646):0.19033808188813469697,1:0.11198373146590309690):0.0;\n";
//
//        trees[3] =    "(12:0.03074164853080323806,(14:0.28354235267611582749,((((6:0.24575236511190479138,3:0.27669877750004706618):0.12346183000473436031,(16:0.27345323581139363567,(2:0.25952584041762544187,(7:0.13314735215857465001,5:0.10564854402522012855):0.10104371710261253470):0.05362992296414054172):0.03200273217217284444):0.36863248264167475243,(((15:0.02110223988704504339,4:0.00639777557208858969):0.24768510157235693048,10:0.25756028527529672978):0.20092006982273108440,((O:2.56103098572207965589,8:0.00000100000050002909):0.16812944217703390537,(9:0.07433198611216682017,13:0.06431981546756693202):0.15576097960789220265):0.23340136310324355318):0.39182426082511839782):0.06886270571453080580,11:0.74081186616078464890):0.49783777059213268901):0.22400130822912850403,1:0.03845212136032683664):0.0;";
//
//        trees[4] =     "((14:0.34535975075053071404,((((8:0.18446096698659020241,(9:0.07574019887156854813,(O:3.02485274708579332525,13:0.00000100000050002909):0.06031186103974060103):0.13204997139393842565):0.14045148750977609620,(10:0.21012017559197199157,(15:0.02693202332471199625,4:0.01312345369941220544):0.19509727373945162654):0.22144481936340332240):0.50440371637957892581,(((16:0.22089459168050687432,(5:0.11546718620580789794,7:0.14793454828801871903):0.04465936833361632691):0.03519308103446531361,2:0.20857735829220303003):0.07125475708349070236,(3:0.24219991474567206802,6:0.22686602857724511018):0.14713876152871011027):0.40482303499946870051):0.10593608537552247351,11:0.70645159295127057941):0.39240867219092495644):0.22086743505460210146,12:0.04368046736286994303,1:0.06694673556289100491):0.0;";
//
//        trees[5] = "((12:0.01244124831499836617,(14:0.36068117687360867762,(((3:0.15050128895012682140,6:0.19737265897379524526):0.19606471557493818092,((16:0.23067097411708073484,2:0.27093012233911950570):0.04796742022486685103,(5:0.12934022296182287604,7:0.11667878826659762448):0.20893031122751259954):0.04748429350467037874):0.39868067894300274201,((((15:0.10740043276400766670,4:0.09203335437680774944):0.12862692486960666804,10:0.23034062124414150796):0.24234543210116527012,(8:0.15275556531350656386,(9:0.03180819295747101122,13:0.03184579283397081378):0.11388367986567096701):0.27587422400224498320):0.45017989548102682607,11:0.77805197595362463048):0.02523327406088173366):0.48066835450698136567):0.35302559292526503043):0.02026037738338782451,O:3.60345000310469165683,1:0.00000100000050002909):0.0;";
//
//        trees[6] =   "(12:0.02665426622803612408,(14:0.26913387902274188423,(11:0.70095360764794667485,((((9:0.08427260969390176770,(O:3.05208532199930759887,13:0.00000100000050002909):0.12191381078469575028):0.06783343721172939433,8:0.18666082839486772271):0.29231797562118899680,(10:0.24336762113110613215,(4:0.07930047814945728746,15:0.07130201072203880275):0.19138898577696586423):0.10804350562318890128):0.61687235837296161556,((6:0.25574059155492367967,3:0.23355128825223403499):0.15823278488642877893,(16:0.25619133815335526361,((7:0.12883286248509762117,5:0.16950675695960301659):0.10915218646184607354,2:0.21050128984494193007):0.02226819916973757363):0.09372695125844936193):0.41291379862869354023):0.06151120265282518373):0.42476731839700593296):0.28998410392102430544,1:0.02928983239912606382):0.0;";
//
//        trees[7] =   "(12:0.01543794591887487182,((11:0.84380373806123631919,((((16:0.22265151229826210533,(7:0.15051539001411451402,5:0.13515038681168714541):0.06413432121655598950):0.01489176281347295676,2:0.22213942552089518956):0.23833638524685848981,(3:0.20713706372306864290,6:0.24561605401480843258):0.12901436615984873679):0.48321801163551925118,((8:0.14423931099104453368,(9:0.09405489340426907707,13:0.10278431425142579847):0.02862701204941410504):0.34843669845758823422,(10:0.23013505224915065117,(15:0.01616787196750350472,4:0.01552249489584072405):0.20316732424115255262):0.19707714640911810910):0.27499432518503919010):0.06759436990636860654):0.45546543422530999523,(O:2.45258603865712565550,14:0.00000100000050002909):0.25942041199689297448):0.35601573847302570019,1:0.03410283571940193187):0.0;";
//
//        trees[8] =    "(12:0.01892630689498123187,(14:0.39816490099406076641,(11:0.73690521736460168345,(((6:0.22722375016873497255,3:0.19401272153340656801):0.09291032803184033628,((2:0.22635049292556133627,(7:0.17822677001188630408,5:0.15535628785717209288):0.06866039891950555252):0.07508973625529276330,16:0.29543954893762708336):0.07915398826744531047):0.42440804379772800736,((8:0.17268201530331817040,((O:3.15399315320461104406,13:0.00000100000050002909):0.04519547417514546001,9:0.08316459542808538052):0.03513243799442695992):0.29138659574922359807,(10:0.30384206006800379418,(15:0.00690278626525097413,4:0.00724117284784580630):0.26414813132350151159):0.13728049462168190331):0.62317981140252709515):0.12219246855830499798):0.35312245113318563305):0.26253761332975883436,1:0.01794494118807227895):0.0;";
//
//        trees[9] =     "(12:0.02204864678559188135,(14:0.30965174234718972901,(((((15:0.01398761381174369427,4:0.00319890353523479389):0.30111727117179659663,(O:2.46592097094110807021,10:0.07643937156956069623):0.22236479283431648257):0.16615258777713340366,(8:0.11784612205464572998,(9:0.11072427820763887374,13:0.12053554257610800404):0.00545785274608394377):0.37270223483332576286):0.42974832440887111851,((3:0.17537983386562600030,6:0.20712498910484364312):0.31904193691708587766,((5:0.11050921647466772912,7:0.11324804800455570697):0.16115111168390866547,(16:0.24096122312784845620,2:0.27449981612747925697):0.01727490839090758262):0.13671665017546061338):0.24242114578412660153):0.14388150264806404133,11:0.64909436999748526542):0.38373014418896267408):0.28473061184493725362,1:0.03619123576604073994):0.0;";


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
//        trees[0] = "(15:0.06385239050018931550,((((7:0.35507858635141326120,(5:0.29353650672204489869,2:0.23406687701105596822):0.06643191841041275192):0.43443454037783985067,((14:0.12499816664032974145,9:0.14170502523795053262):0.04374612299215879102,13:0.15473405018419678081):0.52957691547077978544):0.53786789459448192119,(8:0.12080737028319389614,((12:0.08356589617662767144,11:0.06447326508961689906):0.03147560317215407355,((O:2.34057499021339276979,4:0.00000100000050002909):0.05174336931150425728,16:0.05053598520394290278):0.04770292865750219713):0.03704580009898181259):0.14173872698785303093):0.05301445258794363485,((3:0.06392832764165792825,10:0.05796107858634157867):0.12011281425939168699,6:0.20601493580392435390):0.07731084362708398605):0.24172597823533428785,1:0.08367374764129421882):0.0;";
//
//        trees[1] =   "(((3:0.07067676949193306413,10:0.05490362740030942101):0.15915501529849843121,(6:0.18583036929393462189,((((12:0.03365158635464698955,8:0.05339983396435030866):0.07290109521790824609,11:0.08699423187957987247):0.30145138455352843421,((((O:2.79339209973832547362,5:0.00000100000050002909):0.32374386394702830971,2:0.30392381730752554558):0.04923391482203538783,7:0.38808719386447132482):0.41867808904179854013,((13:0.13693194675898351376,14:0.13127426326912855492):0.05582651902402614025,9:0.20796222415194326461):0.64009955891090419833):0.52676112470034985868):0.10992175270616143246,(16:0.07230998576624945995,4:0.06586836885251330653):0.23735830000761412584):0.05517271419130038296):0.02295560794657667755):0.23447979019970513481,15:0.00599934813918216042,1:0.00406576464609070994):0.0;";
//
//        trees[2] =   "(((6:0.21864692152097200961,(((16:0.09948079467352581584,(12:0.06811430488209768708,8:0.07431607512919138903):0.02648238011667029476):0.03778550672047790993,(4:0.13060908038257934560,11:0.08562746962628120517):0.01219604300873776588):0.14581266225694869410,((13:0.16237848382081532250,(O:2.34077418772927670076,(14:0.14089121826408579352,9:0.11920005236362521406):0.16860184982863507530):0.02427074271366035096):0.59314355288001452404,((5:0.21582170498353936416,2:0.23004171136387729923):0.10360571271432542084,7:0.35052215593934077376):0.61248499378284471106):0.65216579736076796259):0.05585956028786604605):0.00668098040161886706,(3:0.14093840019188050294,10:0.14105455128389468578):0.06270524915762901164):0.15780604191116268109,15:0.04601464874506398101,1:0.03925786701673324858):0.0;";
//
//        trees[3] =    "(15:0.01504073197183840652,((((((2:0.27715353055650876479,5:0.23994623184437588459):0.16960311258998658768,7:0.27429391383418738615):0.57488257589166602557,(13:0.09258402465310638929,(14:0.09665510942607179190,(O:3.03517750932576335643,9:0.00000100000050002909):0.11404251417786588629):0.09216691812612336221):0.56963286005863877914):0.43066569270239335454,(11:0.13523535475148720786,(12:0.09422117869190024808,((16:0.06566091067862804553,4:0.06232752953366831050):0.02686209599698902856,8:0.08478586141789691555):0.01379996860303076521):0.03761460486523916097):0.17379066189132019726):0.10499904924467018397,(3:0.04023156769678863653,10:0.03203571117735003887):0.23644419357652782243):0.01263779425747846991,6:0.22461533901501193333):0.21618297005840378389,1:0.01669062174999247361):0.0;";
//
//        trees[4] =     "((((((9:0.14346674889858623825,14:0.15915614783964654455):0.05233336282342995749,13:0.17119692549367435119):0.45529474007969017446,((2:0.15830736928244024120,5:0.19584476185184726549):0.25736835478331748783,(O:2.48142909549519119494,7:0.00000100000050002909):0.36180239062652264082):0.61580560652742155625):0.52341327161783202104,(11:0.09770273434829912507,(12:0.09415447162603486631,(8:0.08649011051468645517,(4:0.06468901630280333992,16:0.07160461603597464975):0.02931889857258338991):0.00345164581426198741):0.01330590757513378138):0.22187787603728831298):0.12127200951536722118,((10:0.04019493585650776163,3:0.04866162700196443452):0.17646419724619608060,6:0.21046112688285797954):0.01347532607852260957):0.21839907362733285145,15:0.02087511746025366710,1:0.02122372485415865984):0.0;";

//        trees[5] = "(15:0.06385239050018931550,((((7:0.35507858635141326120,(5:0.29353650672204489869,2:0.23406687701105596822):0.06643191841041275192):0.43443454037783985067,((14:0.12499816664032974145,9:0.14170502523795053262):0.04374612299215879102,13:0.15473405018419678081):0.52957691547077978544):0.53786789459448192119,(8:0.12080737028319389614,((12:0.08356589617662767144,11:0.06447326508961689906):0.03147560317215407355,((O:2.34057499021339276979,4:0.00000100000050002909):0.05174336931150425728,16:0.05053598520394290278):0.04770292865750219713):0.03704580009898181259):0.14173872698785303093):0.05301445258794363485,((3:0.06392832764165792825,10:0.05796107858634157867):0.12011281425939168699,6:0.20601493580392435390):0.07731084362708398605):0.24172597823533428785,1:0.08367374764129421882):0.0;";

//        trees[6] =   "(((3:0.07067676949193306413,10:0.05490362740030942101):0.15915501529849843121,(6:0.18583036929393462189,((((12:0.03365158635464698955,8:0.05339983396435030866):0.07290109521790824609,11:0.08699423187957987247):0.30145138455352843421,((((O:2.79339209973832547362,5:0.00000100000050002909):0.32374386394702830971,2:0.30392381730752554558):0.04923391482203538783,7:0.38808719386447132482):0.41867808904179854013,((13:0.13693194675898351376,14:0.13127426326912855492):0.05582651902402614025,9:0.20796222415194326461):0.64009955891090419833):0.52676112470034985868):0.10992175270616143246,(16:0.07230998576624945995,4:0.06586836885251330653):0.23735830000761412584):0.05517271419130038296):0.02295560794657667755):0.23447979019970513481,15:0.00599934813918216042,1:0.00406576464609070994):0.0;";
//
//        trees[7] =   "(((6:0.21864692152097200961,(((16:0.09948079467352581584,(12:0.06811430488209768708,8:0.07431607512919138903):0.02648238011667029476):0.03778550672047790993,(4:0.13060908038257934560,11:0.08562746962628120517):0.01219604300873776588):0.14581266225694869410,((13:0.16237848382081532250,(O:2.34077418772927670076,(14:0.14089121826408579352,9:0.11920005236362521406):0.16860184982863507530):0.02427074271366035096):0.59314355288001452404,((5:0.21582170498353936416,2:0.23004171136387729923):0.10360571271432542084,7:0.35052215593934077376):0.61248499378284471106):0.65216579736076796259):0.05585956028786604605):0.00668098040161886706,(3:0.14093840019188050294,10:0.14105455128389468578):0.06270524915762901164):0.15780604191116268109,15:0.04601464874506398101,1:0.03925786701673324858):0.0;";
//
//        trees[8] =    "(15:0.01504073197183840652,((((((2:0.27715353055650876479,5:0.23994623184437588459):0.16960311258998658768,7:0.27429391383418738615):0.57488257589166602557,(13:0.09258402465310638929,(14:0.09665510942607179190,(O:3.03517750932576335643,9:0.00000100000050002909):0.11404251417786588629):0.09216691812612336221):0.56963286005863877914):0.43066569270239335454,(11:0.13523535475148720786,(12:0.09422117869190024808,((16:0.06566091067862804553,4:0.06232752953366831050):0.02686209599698902856,8:0.08478586141789691555):0.01379996860303076521):0.03761460486523916097):0.17379066189132019726):0.10499904924467018397,(3:0.04023156769678863653,10:0.03203571117735003887):0.23644419357652782243):0.01263779425747846991,6:0.22461533901501193333):0.21618297005840378389,1:0.01669062174999247361):0.0;";
//
//        trees[9] =     "((((((9:0.14346674889858623825,14:0.15915614783964654455):0.05233336282342995749,13:0.17119692549367435119):0.45529474007969017446,((2:0.15830736928244024120,5:0.19584476185184726549):0.25736835478331748783,(O:2.48142909549519119494,7:0.00000100000050002909):0.36180239062652264082):0.61580560652742155625):0.52341327161783202104,(11:0.09770273434829912507,(12:0.09415447162603486631,(8:0.08649011051468645517,(4:0.06468901630280333992,16:0.07160461603597464975):0.02931889857258338991):0.00345164581426198741):0.01330590757513378138):0.22187787603728831298):0.12127200951536722118,((10:0.04019493585650776163,3:0.04866162700196443452):0.17646419724619608060,6:0.21046112688285797954):0.01347532607852260957):0.21839907362733285145,15:0.02087511746025366710,1:0.02122372485415865984):0.0;";

//

//        trees[0] = "((((11:0.70644398931864882396,(3:0.34767924403943067935,(((5:0.15710353984118755055,7:0.18935903425796674071):0.08874855967111235089,(16:0.20087264366883630928,2:0.25783647194418191084):0.03843178499567419248):0.05896998829160856465,6:0.38085907307506022823):0.07311579398713173583):0.37200568204787165527):0.01921612205860232306,(((13:0.07504554726485214433,9:0.06313470059355173747):0.06576799914386330381,8:0.14254521883365151025):0.26792964654666351354,(10:0.25828580936016115599,(15:0.01424732461309747601,4:0.01322523820803084803):0.18785549914016130835):0.19246845125134842691):0.55215562996529365947):0.46131833154437729139,(O:4.05911231085117574224,14:0.00000100000050002909):0.29596406113180256936):0.28137585062420888971,12:0.01731115410631995505,1:0.02501476357694640967):0.0;";
//
//        trees[1] =   "(12:0.04399827182659500624,(((((10:0.25444716531487121536,(15:0.00355401690007028448,4:0.00549988354214483999):0.19366689258000677198):0.18235892591196123380,(8:0.11046727708066554630,(9:0.07488860481588313633,13:0.04400698866155940830):0.07220305126224457748):0.22117053026417671235):0.54084844484732019776,((2:0.25676556371336750617,(16:0.19503163588057048639,(5:0.14646252764275405345,7:0.19104839503609721829):0.08411985892033234380):0.06458311226381635572):0.07172982922310700005,(3:0.30548101010256234478,6:0.31611224179943486767):0.03510270822114055361):0.36270266518108468246):0.00000100000050002909,(O:2.77115580175478282143,11:0.00000100000050002909):0.62484388145985103602):0.43031475524696627399,14:0.25858177305660351442):0.27488441331466606510,1:0.03821119722230179039):0.0;";
//
//        trees[2] =   "(12:0.05389264015061668439,(14:0.31036431233883499115,((((8:0.12617679588613076036,(13:0.05779635336035422732,(O:2.36442923010710437381,9:0.00000100000050002909):0.05984303062817749458):0.07391840543060251845):0.34912355218082097297,(10:0.20916395904115100968,(4:0.02087354149109655471,15:0.03389591261305226894):0.16063257781127360246):0.19088129202874004609):0.48533178418738748272,((2:0.19709442058381326790,(16:0.22116446851592569400,(7:0.14106565037434873333,5:0.14933792283054186933):0.09067034675596909177):0.02074605264692864354):0.20315483363510919435,(6:0.17929711217024171810,3:0.22547757042749946854):0.03814184753746297646):0.54190226689483522726):0.06556032698135033643,11:0.74641835382036614277):0.34983134815328109646):0.19033808188813469697,1:0.11198373146590309690):0.0;\n";
//
//        trees[3] =    "(12:0.03074164853080323806,(14:0.28354235267611582749,((((6:0.24575236511190479138,3:0.27669877750004706618):0.12346183000473436031,(16:0.27345323581139363567,(2:0.25952584041762544187,(7:0.13314735215857465001,5:0.10564854402522012855):0.10104371710261253470):0.05362992296414054172):0.03200273217217284444):0.36863248264167475243,(((15:0.02110223988704504339,4:0.00639777557208858969):0.24768510157235693048,10:0.25756028527529672978):0.20092006982273108440,((O:2.56103098572207965589,8:0.00000100000050002909):0.16812944217703390537,(9:0.07433198611216682017,13:0.06431981546756693202):0.15576097960789220265):0.23340136310324355318):0.39182426082511839782):0.06886270571453080580,11:0.74081186616078464890):0.49783777059213268901):0.22400130822912850403,1:0.03845212136032683664):0.0;";
//
//        trees[4] =     "((14:0.34535975075053071404,((((8:0.18446096698659020241,(9:0.07574019887156854813,(O:3.02485274708579332525,13:0.00000100000050002909):0.06031186103974060103):0.13204997139393842565):0.14045148750977609620,(10:0.21012017559197199157,(15:0.02693202332471199625,4:0.01312345369941220544):0.19509727373945162654):0.22144481936340332240):0.50440371637957892581,(((16:0.22089459168050687432,(5:0.11546718620580789794,7:0.14793454828801871903):0.04465936833361632691):0.03519308103446531361,2:0.20857735829220303003):0.07125475708349070236,(3:0.24219991474567206802,6:0.22686602857724511018):0.14713876152871011027):0.40482303499946870051):0.10593608537552247351,11:0.70645159295127057941):0.39240867219092495644):0.22086743505460210146,12:0.04368046736286994303,1:0.06694673556289100491):0.0;";
//
//        trees[5] = "((12:0.01244124831499836617,(14:0.36068117687360867762,(((3:0.15050128895012682140,6:0.19737265897379524526):0.19606471557493818092,((16:0.23067097411708073484,2:0.27093012233911950570):0.04796742022486685103,(5:0.12934022296182287604,7:0.11667878826659762448):0.20893031122751259954):0.04748429350467037874):0.39868067894300274201,((((15:0.10740043276400766670,4:0.09203335437680774944):0.12862692486960666804,10:0.23034062124414150796):0.24234543210116527012,(8:0.15275556531350656386,(9:0.03180819295747101122,13:0.03184579283397081378):0.11388367986567096701):0.27587422400224498320):0.45017989548102682607,11:0.77805197595362463048):0.02523327406088173366):0.48066835450698136567):0.35302559292526503043):0.02026037738338782451,O:3.60345000310469165683,1:0.00000100000050002909):0.0;";
//
//        trees[6] =   "(12:0.02665426622803612408,(14:0.26913387902274188423,(11:0.70095360764794667485,((((9:0.08427260969390176770,(O:3.05208532199930759887,13:0.00000100000050002909):0.12191381078469575028):0.06783343721172939433,8:0.18666082839486772271):0.29231797562118899680,(10:0.24336762113110613215,(4:0.07930047814945728746,15:0.07130201072203880275):0.19138898577696586423):0.10804350562318890128):0.61687235837296161556,((6:0.25574059155492367967,3:0.23355128825223403499):0.15823278488642877893,(16:0.25619133815335526361,((7:0.12883286248509762117,5:0.16950675695960301659):0.10915218646184607354,2:0.21050128984494193007):0.02226819916973757363):0.09372695125844936193):0.41291379862869354023):0.06151120265282518373):0.42476731839700593296):0.28998410392102430544,1:0.02928983239912606382):0.0;";
//
//        trees[7] =   "(12:0.01543794591887487182,((11:0.84380373806123631919,((((16:0.22265151229826210533,(7:0.15051539001411451402,5:0.13515038681168714541):0.06413432121655598950):0.01489176281347295676,2:0.22213942552089518956):0.23833638524685848981,(3:0.20713706372306864290,6:0.24561605401480843258):0.12901436615984873679):0.48321801163551925118,((8:0.14423931099104453368,(9:0.09405489340426907707,13:0.10278431425142579847):0.02862701204941410504):0.34843669845758823422,(10:0.23013505224915065117,(15:0.01616787196750350472,4:0.01552249489584072405):0.20316732424115255262):0.19707714640911810910):0.27499432518503919010):0.06759436990636860654):0.45546543422530999523,(O:2.45258603865712565550,14:0.00000100000050002909):0.25942041199689297448):0.35601573847302570019,1:0.03410283571940193187):0.0;";
//
//        trees[8] =    "(12:0.01892630689498123187,(14:0.39816490099406076641,(11:0.73690521736460168345,(((6:0.22722375016873497255,3:0.19401272153340656801):0.09291032803184033628,((2:0.22635049292556133627,(7:0.17822677001188630408,5:0.15535628785717209288):0.06866039891950555252):0.07508973625529276330,16:0.29543954893762708336):0.07915398826744531047):0.42440804379772800736,((8:0.17268201530331817040,((O:3.15399315320461104406,13:0.00000100000050002909):0.04519547417514546001,9:0.08316459542808538052):0.03513243799442695992):0.29138659574922359807,(10:0.30384206006800379418,(15:0.00690278626525097413,4:0.00724117284784580630):0.26414813132350151159):0.13728049462168190331):0.62317981140252709515):0.12219246855830499798):0.35312245113318563305):0.26253761332975883436,1:0.01794494118807227895):0.0;";
//
//        trees[9] =     "(12:0.02204864678559188135,(14:0.30965174234718972901,(((((15:0.01398761381174369427,4:0.00319890353523479389):0.30111727117179659663,(O:2.46592097094110807021,10:0.07643937156956069623):0.22236479283431648257):0.16615258777713340366,(8:0.11784612205464572998,(9:0.11072427820763887374,13:0.12053554257610800404):0.00545785274608394377):0.37270223483332576286):0.42974832440887111851,((3:0.17537983386562600030,6:0.20712498910484364312):0.31904193691708587766,((5:0.11050921647466772912,7:0.11324804800455570697):0.16115111168390866547,(16:0.24096122312784845620,2:0.27449981612747925697):0.01727490839090758262):0.13671665017546061338):0.24242114578412660153):0.14388150264806404133,11:0.64909436999748526542):0.38373014418896267408):0.28473061184493725362,1:0.03619123576604073994):0.0;";
        List tempGTs = new ArrayList<Tree>();
        double ogMin = Double.POSITIVE_INFINITY;
        double ogMax = 0.0;
        for(int i=0;i<lociNum;i++){
            Tree tempT = Trees.readTree(trees[i]);
            System.out.println("# " + i + operator.getDistance(tempT, Trees.readTree(trueGTS.get(i))));
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
        //        currentHGTLL[0] = -16449.806508;
//        currentHGTLL[1] = -16173.220201;
//        currentHGTLL[2] = -16446.524409;
//        currentHGTLL[3] = -16562.363070;
//        currentHGTLL[4] = -16412.064862;
//        currentHGTLL[5] = -16696.373035;
//        currentHGTLL[6] = -17075.786169;
//        currentHGTLL[7] = -16473.223523;
//        currentHGTLL[8] = -16395.286595;
//        currentHGTLL[9] = -16610.970346;
//        currentHGTLL[5] = -14669.564935;
//        currentHGTLL[6] = -14539.299872;
//        currentHGTLL[7] = -14887.092847;
//        currentHGTLL[8] = -14644.965958;
//        currentHGTLL[9] = -15246.248191;

//        currentHGTLL[0] = -14885.327181;
//        currentHGTLL[1] = -14376.336898;
//        currentHGTLL[2] = -15231.805277;
//        currentHGTLL[3] = -14575.994606;
//        currentHGTLL[4] = -14400.796516;
//        currentHGTLL[5] = -14776.970112;
//        currentHGTLL[6] = -14702.383736;
//        currentHGTLL[7] = -14283.969575;
//        currentHGTLL[8] = -14642.651198;
//        currentHGTLL[9] = -14987.974635;
//        currentHGTLL[0] = -16449.806508;
//        currentHGTLL[1] = -16173.220201;
//        currentHGTLL[2] = -16446.524409;
//        currentHGTLL[3] = -16562.363070;
//        currentHGTLL[4] = -16412.064862;
//        currentHGTLL[5] = -16696.373035;
//        currentHGTLL[6] = -17075.786169;
//        currentHGTLL[7] = -16473.223523;
//        currentHGTLL[8] = -16395.286595;
//        currentHGTLL[9] = -16610.970346;
        //List<Tree> gtsOG = operator.getOGGTS(gts, ogHeight[0]);
        INIT_ST = inferSTByAST(tempGTs);
        //optimizeBL(tempGTs,INIT_ST,ogHeight[0]);
        List<Double> proGT = getGTSLLBySTYF(finalGTS, INIT_ST);
        currentHLL = 0.0;
        currentLL = 0.0;

        for(int i = 0;i<proGT.size();i++) {
            currentHLL += currentHGTLL[i];
            //for init, ProGT = -INFINITY
            if(proGT.get(i)== Double.NEGATIVE_INFINITY)
                proGT.set(i,currentHGTLL[i]);
            currentGTLL[i] = currentHGTLL[i] + proGT.get(i);
            currentLL += currentGTLL[i];
        }
        if(FULL_LL)
            GLOBALBESTLL = currentLL;
        else
            GLOBALBESTLL = currentHLL;

        double p2 = currentHLL-currentLL;
        if (FULL_LL)
            iter_LL = currentHLL*ll_Ratio + p2*(2-ll_Ratio);
        else
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
            msDist_Locus[0][k] = 1;
            thisD += d;
        }
        gtDis[0] = thisD / lociNum;
        return INIT_ST;
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

//    public List<List<Tree>> updateOPTGTS(List<String> topos, List<String> cTopos, List<List<Tree>> optGTList) throws IOException {
//        List<List<Tree>> optGTSi = new ArrayList<List<Tree>>();
//        for(int i =0;i<topos.size();i++) {
//            String temp = topos.get(i);
//            if(cTopos.contains(temp)){
//                int j = 0;
//                for(;j<cTopos.size();j++){
//                    if(cTopos.get(i).equals(temp))
//                        break;
//                }
//                optGTSi.add(optGTList.get(j));
//                continue;
//            }
//
//            else {
//                List<Tree> tempList = operator.optimizeBLGT(temp, trueSeq, halfTheta);
//                optGTList.add(tempList);
//                cTopos.add(temp);
//                optGTSi.add(tempList);
//            }
//        }
//        return optGTSi;
//    }


//    public Double getBestGTS(List<String> topos, Network tempST, List<Tree> bestGTSi ) throws IOException {
//        List<List<Tree>> optGTSi = updateOPTGTS(topos, checkedTopo, optGTSList);
//        Double llForIter_i = 0.0;
//        for(int i = 0;i<lociNum;i++){
//            Double iLL = getBestGTi(optGTSi.get(i), bestGTSi, tempST);
//            currentGTLL[i] = iLL;
//            llForIter_i+=iLL;
//        }
//        return llForIter_i;
//    }


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

//    public List<List<Tree>> updateOPTGTS(List<String> topos, List<String> cTopos, List<List<Tree>> optGTList) throws IOException {
//        List<List<Tree>> optGTSi = new ArrayList<List<Tree>>();
//        for(int i =0;i<topos.size();i++) {
//            String temp = topos.get(i);
//            List<Tree> tempList = operator.optimizeBLGT(temp, trueSeq, halfTheta);
//            //optGTList.add(tempList);
//            cTopos.add(temp);
//            optGTSi.add(tempList);
//
//        }
//        return optGTSi;
//    }

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