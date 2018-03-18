package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

/**
 * Created by doriswang on 3/18/18.
 * for bio dataset
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
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.TemporalConstraints;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GLASSInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.*;
public class BioInfer48 {

    private static String OutDir;
    protected static String _msdir;
    protected static String _ASTRALdir;
    private static String RESULT_DIR;
    private int TREEINDEX;

    private Network<NetNodeInfo> INIT_ST; //cUnit
    private Network<NetNodeInfo> currentST;
    private static String trueST;
    private Tree Constraint_ST;// GLASS tree from UPGMA. For external branch

    private List<String> trueGTS;
    private List<Tree> trueTopos;
    private List<Tree> finalGTS;
    private static List<Alignment> trueSeq;


    private String GLOBALBESTST;
    private List<String> GLOBALBESTGTS;


//    //7.26
//    private List<String> checkedTopo;  //topos as checked order
//    private List<List<Tree>> optGTSList;  // For one locus l, locus<map<bestGT,LL>>
//

    private IIGTSimulator simulator;
    private InferOperator operator;

    //TODO: PARA
    private static int ITERATION; //{1k, 5k ,10k}
    private static boolean ISREALDATA;
    private static boolean FULL_LL;
    private static boolean GET_MS_DETAIL;

    private static double halfTheta;//need to init in ALL classes ||  theta -> POP_SIZE
    private static Integer lociNum;
    private static Integer refineGTSize;

    private static double[] _scales; //temp never used
    private static int[] _seqLens;
    private static int taxaNum;
    private int itNum;// iteration number
    private int endCounter;// iteration number
    private List<Tree> locusBestGTS;
    private double[] weights;
    private double ll_Ratio; // (For Expectation) : E = p(Seq|raxml_G) * ratio P(ms_G|ST)   because P(G|ST) is not accurate
    private boolean ll_Ratio_Change;
    private static String STRATEGY; //RANDOM IMPROVE NONE
    private List<Integer> bestITNum;///temporary best iteration number for iteration i


    //likelihood:
    private double GLOBALBESTLL;
    private double GLOBALBESTHLL;// ll(Seq|GT)
    private double GLOBALBESTSTLL;// ll(GT|ST)
    private double currentLL;//=gtsLL[iter] LL(Seq|GTS)*(GTS|ST) for all locus for iteration i
    private double currentHLL;// =gtsHLL[iter]
    private double currentSTLL;// = current ll(G|ST)
    private int bestLLIter;//iteration number for best Global ST
    private double maxLL; //current best LL

    private double[] locusGTLL; // P(Seq|GT)
    private double[] iterLLSeq;// LL(Seq|GTS) for each locus for iteration i
    private double[] iterLLST;// LL(Seq|GTS) for each locus for iteration i
    private double[] gtsLL;// LL(Seq|GTS)*(GTS|ST) for all locus for iteration i
    private double[] gtsHLL;// LL(Seq|GTS) for all locus for iteration i
    private double[] currentITLL; // tempLL for iteration i (may be rejected)
    private double[] currentITHLL; //temp LL(Seq|GT) for iteration i (may be rejected)
    private double[] currentITSTLL; //temp LL(GT|ST) for iteration i (may be rejected)
    private List<Double> bestSTLL;// temporary best ST LL for iteration i

    //double bestLL = llList[0] + ll_Ratio*msTopoLL.get(i);

    //distance:
    private  double[] astInit; //ASTRAL : [0]unrooted_ST dist [1]unrooted_GT dist
    private double[][] gtDis_Locus;// gt distance for each locus for each iteration
    private double[][] msDist_Locus;// ms_gt distance for each locus for each iteration --> to test the efficiency of MS
    private List<Double> bestSTD;//temporary  ST distance for GB_bestST for iteration i
    private double[] stDis;//ST_i for each iteration
    private double[] stRootedDis;//ST_i for each iteration
    private double[] gtDis; // GT_i for temp result(include rej)
//

    //All parameters should be initialed here
    BioInfer48(String inputPath, String outputPath, int index, double hTheta, int iteration, int lnum, boolean msDetail, int rNum, int msSize) throws IOException, InterruptedException, ParseException {

        String fTemp = "syr";
        TREEINDEX = index*lnum;
        endCounter = 0;
        GET_MS_DETAIL = msDetail;
        //index for program   TREEINDEX for inputData
        _msdir = inputPath + "tools/msFiles/msdir/ms"; ///home/yw58/IIG/tools/Astral/astral.4.11.1.jar
        _ASTRALdir = inputPath + "tools/Astral/" + rNum+ "/";
        RESULT_DIR = outputPath ;//"/Users/doriswang/PhyloNet/Data/IIG/result/" | home/yw58/IIG/output
        trueGTS = new ArrayList<String>();
        finalGTS = new ArrayList<Tree>();
        trueTopos = new ArrayList<Tree>();


        GLOBALBESTGTS = new ArrayList<String>();
        trueSeq = new ArrayList<Alignment>();
        //TODO: PARA
        ITERATION = iteration; //{1k, 5k ,10k}
        ISREALDATA = false;
        FULL_LL = true;
        halfTheta = hTheta;//need to init in ALL classes ||  theta -> POP_SIZE
        lociNum = lnum;
        //ifOutGroup = false;

        //ogHeight = new double[3]; // og height from all initGT[min, max]
        refineGTSize = msSize;
        _scales = new double[]{1.0};
        _seqLens = new int[]{471,525,420,260,260,469,796,296,220,262,256,194,447,849,798,522,267,274,665};
        weights = new double[lnum];
        for(int i =0;i<lnum;i++){
            weights[i] = 1;
        }

        taxaNum = 48;
        ll_Ratio = 1.0; // P(seq|GT): p(GT|ST)
        ll_Ratio_Change = false;
        gtsLL = new double[ITERATION];// LL(Seq|GTS)*(GTS|ST)
        // LL(Seq|GTS)*(GTS|ST) for pseudo version || not used for classic version
        gtsHLL = new double[ITERATION];
        currentITLL = new double[ITERATION];
        currentLL = 0.0; // LL(Seq|GTS)*(GTS|ST)
        currentHLL = 0.0;
        currentITHLL = new double[ITERATION];
        currentITSTLL = new double[ITERATION];
        gtDis = new double[ITERATION];
        stDis = new double[ITERATION];
        stRootedDis = new double[ITERATION];
        gtDis_Locus = new double[ITERATION][lociNum];
        msDist_Locus = new double[ITERATION][lociNum];
        bestLLIter = 0;
        bestITNum = new ArrayList<Integer>();
        bestSTLL = new ArrayList<Double>();
        bestSTD = new ArrayList<Double>();
        STRATEGY = "RANDOM"; //RANDOM IMPROVE N
        //TODO Simulation
        simulator = new IIGTSimulator(lociNum, _scales, _seqLens, halfTheta, ITERATION, inputPath, rNum);
        operator = new InferOperator(inputPath,outputPath,index, rNum);
        operator.changeUpdateSh(lociNum,refineGTSize);
        OutDir = RESULT_DIR  + fTemp + "/" + lociNum + "/" + ITERATION +"/"+ halfTheta*100000 + "/" + _seqLens[0] + "/" + msSize + "/"; // output/0001/1/
        operator.isExitsPath(OutDir);
        locusBestGTS = new ArrayList<Tree>();
        locusGTLL = new double[lociNum];
        //{1,2,3,4,5,1,2,3,4,5};
        //{1,1,1,1,1,5,5,5,5,5};
        astInit = new double[6];
        //String tempPath = inputPath + "input/" + fTemp + "/38LocusData/";
        //"a/input/"
        if(fTemp.equals("syr")){
            trueST = "(((((24,50),(25,51)),(20,46)),((48,((22,23),49)),(21,47))),(((((36,(10,(7,33))),((3,29),(4,30))),((2,28),(6,32))),(((5,31),(8,34)),(9,35))),((((((13,39),(14,40)),(11,37)),(12,38)),(((16,42),(17,43)),(15,41))),((18,44),(19,45)))));";
            String tempPath = inputPath + "input/" + fTemp + "/38LocusData/";
            trueSeq = operator.loadSYRSeq48(lociNum, tempPath, taxaNum, trueSeq);
            //TODO 9-11 load startTree?
            for(int i = 0;i<lociNum;i++){
                trueGTS.add("(((((24,50),(25,51)),(20,46)),((48,((22,23),49)),(21,47))),(((((36,(10,(7,33))),((3,29),(4,30))),((2,28),(6,32))),(((5,31),(8,34)),(9,35))),((((((13,39),(14,40)),(11,37)),(12,38)),(((16,42),(17,43)),(15,41))),((18,44),(19,45)))));");
            }
        }

        for(int i = 0; i<lociNum; i++){
            locusGTLL[i] = Double.NEGATIVE_INFINITY;
        }
        //weights = {5,1,1,1,1};
    }

    public static void main(String[] args) throws IOException, ParseException, InterruptedException {
//UGARInfer(String inputPath, String outputPath, int index, double hTheta, int iteration, int lnum, boolean msDetail, int rNum, int seqLength, int msSize)
        //BioInfer48 ui = new BioInfer48(args[0],args[1],Integer.valueOf(args[2]),Double.valueOf(args[3]),Integer.valueOf(args[4]), Integer.valueOf(args[5]),Boolean.valueOf(args[6]),Integer.valueOf(args[7]),Integer.valueOf(args[8]));
        BioInfer48 ui = new BioInfer48("/Users/doriswang/PhyloNet/","/Users/doriswang/PhyloNet/Data/IIG/result/8/",0,0.00075,2,19, true , 0,  5);

        // lNum: lociNum:  snack 19
        ui.infer(trueSeq, ITERATION);
        System.out.println("Finish");
    }

    //main function
    //Input: true alignments
    //Output: best ST
    public List<Tree> infer(List<Alignment> aln, int t) throws IOException, ParseException, InterruptedException {
        //CHANGE: PARA
        itNum = 0;
        ITERATION = t;
        long start0 = System.currentTimeMillis();
        List<Tree> gts = operator.initFasttree(trueSeq,_seqLens,lociNum);
        double dist = 0.0;
        for(int i=0;i<gts.size();i++){
            dist+=operator.getDistance(gts.get(i),Trees.readTree(trueGTS.get(i)));
        }
        astInit[1] = dist/gts.size();
        System.out.println("AST  GT distance: " + astInit[1]);

        String n = initAST(gts);
        astInit[0] = operator.getDistance(Trees.readTree(n),Trees.readTree(trueST));
        System.out.println("AST  ST distance: " + astInit[0]);
        System.out.println("AST is " + n);
        long start = System.currentTimeMillis();


        //Part 1: Initialization
        Network<NetNodeInfo> tempST = initST(aln);
        itNum++;
        System.out.println("\n" + "#0" + " iteration start! ");
        System.out.println("Current global best LL is" + GLOBALBESTLL + "  Current global best LL(S|G) is" + GLOBALBESTHLL);
        System.out.println("RF of Best ST : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)));
        System.out.println("Average distance of GTS is " + gtDis[0]);


        //Part 2: Expectation Maximization:
        for (int i = 1; i < t; i++) {
            System.out.println("\n" + "#" + itNum + " iteration start! ");
            //System.out.println(" Current ST is " + tempST.toString());

            //TODO tuning ll ratio
//            if(ll_Ratio_Change) {
//                if (i % 15 == 0) {
//                    ll_Ratio *= 2;
//                    System.out.println("!!!!! ll_Ratio: " + ll_Ratio);
//                }
//            }
            tempST = refineGeneSet(tempST);
            itNum++;
            endCounter++;
            if(endCounter>10)
                break;
        }

        long end = System.currentTimeMillis();
        long costtime = end - start;
        long astTime = start-start0;
        String resultFolder = OutDir ;
        operator.isExitsPath(resultFolder);
        BufferedWriter llOut1 = new BufferedWriter(new FileWriter(resultFolder + "RunningTime.txt"));
        llOut1.write("Running time:" + String.valueOf(costtime) + "\n");
        llOut1.write("AST Running time:" + String.valueOf(astTime) + "\n");
        // System.out.println("------Running time: " + costtime);
        operator.isExitsPath(resultFolder);
        llOut1.flush();
        llOut1.close();
        BufferedWriter llOut = new BufferedWriter(new FileWriter(resultFolder + "acc_LL_Seq.txt"));
        for (int k = 0; k < ITERATION; k++) {
            llOut.write(k + " " + String.valueOf(gtsHLL[k]) + "\n");
        }
        llOut.flush();
        llOut.close();
        BufferedWriter llPOut = new BufferedWriter(new FileWriter(resultFolder + "acc_LL_All.txt"));
        gtsLL[0] = gtsLL[1];
        for (int k = 0; k < ITERATION; k++) {
            llPOut.write(k + " " + String.valueOf(gtsLL[k]) + "\n");
        }
        llPOut.flush();
        llPOut.close();
        BufferedWriter als = new BufferedWriter(new FileWriter(resultFolder + "acc_LL_ST.txt"));

        for (int k = 0; k < ITERATION; k++) {
            als.write(k + " " + String.valueOf(gtsLL[k] - gtsHLL[k]) + "\n");
        }
        als.flush();
        als.close();
        BufferedWriter tlsSeq = new BufferedWriter(new FileWriter(resultFolder + "temp_LL_Seq.txt"));
        currentITHLL[0] = currentITHLL[1];
        for (int k = 0; k < ITERATION; k++) {
            tlsSeq.write(k + " " + String.valueOf(currentITHLL[k]) + "\n");
        }
        tlsSeq.flush();
        tlsSeq.close();
        BufferedWriter tla = new BufferedWriter(new FileWriter(resultFolder + "temp_LL_All.txt"));
        currentITLL[0] = currentITLL[1];
        for (int k = 0; k < ITERATION; k++) {
            tla.write(k + " " + String.valueOf(currentITLL[k]) + "\n");
        }
        tla.flush();
        tla.close();
        BufferedWriter tls = new BufferedWriter(new FileWriter(resultFolder + "temp_LL_ST.txt"));
        currentITSTLL[0] = currentITLL[1];
        for (int k = 0; k < ITERATION; k++) {
            tls.write(k + " " + String.valueOf(currentITLL[k]) + "\n");
        }
        tls.flush();
        tls.close();
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
        for (int tid = 0; tid < lociNum; tid++) {
            String fileName = resultFolder + "GTDistance_" + tid + ".txt";
            BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
            for (int it = 0; it < ITERATION; it++) {
                out.write(it + " " + gtDis_Locus[it][tid] + "\n");
                out.flush();
            }
            out.close();
            if(GET_MS_DETAIL) {
                String fileName1 = resultFolder + "MSGTDistance_" + tid + ".txt";
                BufferedWriter out1 = new BufferedWriter(new FileWriter(fileName1));
                for (int it = 0; it < ITERATION; it++) {
                    out1.write(it + " " + msDist_Locus[it][tid] + "\n");
                    out1.flush();
                }
                out1.close();
            }
        }


        // experiment setting
        BufferedWriter bw2 = new BufferedWriter(new FileWriter(resultFolder + "setting" + ".txt"));
        bw2.write(GLOBALBESTST + "\n");
        bw2.write("STATEGY : " + STRATEGY + "\n");
        bw2.write("Best Likelihood : LL = " + String.valueOf(GLOBALBESTLL) + "  ; P(Seq|GT)" +  String.valueOf(GLOBALBESTHLL) +  "\n");
        bw2.write("Iteration number : " + String.valueOf(bestLLIter) + "\n");
        bw2.write("Likelihood combined ratio : " + ll_Ratio + "\n");
        bw2.write("Number of loci : " + lociNum + "\n");
        bw2.write("Number of refine gene tree set :  " + refineGTSize + "\n");
        bw2.write("Number of iteration :  " + ITERATION + "\n");
        bw2.write("Number of ms Size :  " + refineGTSize + "\n");
        bw2.write("Seq lens are :  ");
        for(int i = 0 ; i< _seqLens.length; i++){
            bw2.write(_seqLens[i] + " ");
        }
        bw2.write("HalfTheta  " + halfTheta + "\n");
        bw2.flush();
        bw2.close();

        // result trees
        String fileName = resultFolder + "GlobalBestTrees.txt";
        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
        out.write("ST: " + GLOBALBESTST + "\n" + "GT: ");
        for (int tid = 0; tid < GLOBALBESTGTS.size(); tid++) {
            String gt = GLOBALBESTGTS.get(tid);
            out.write(gt + "\n");
            out.flush();
        }
        out.close();

        // RAxML trees
        fileName = resultFolder + "RAxMLBestTrees.txt";
        out = new BufferedWriter(new FileWriter(fileName));
        for (int tid = 0; tid < locusBestGTS.size(); tid++) {
            String gt = locusBestGTS.get(tid).toString();
            out.write(gt + "\n");
            out.flush();
        }
        out.close();


        //Distance:
        BufferedWriter gtD = new BufferedWriter(new FileWriter(resultFolder + "GTdistance.txt"));
        for (int gtn = 0; gtn < ITERATION; gtn++) {
            gtD.write(String.valueOf(gtn) + ":" + String.valueOf(gtDis[gtn]) + "\n");
        }
        gtD.flush();
        gtD.close();
        BufferedWriter stD = new BufferedWriter(new FileWriter(resultFolder + "STdistance.txt"));
        for (int stn = 0; stn < ITERATION; stn++) {
            stD.write(String.valueOf(stn) + ":" + String.valueOf(stDis[stn]) + "\n");
        }
        stD.flush();
        stD = new BufferedWriter(new FileWriter(resultFolder + "STRootedDistance.txt"));
        for (int stn = 0; stn < ITERATION; stn++) {
            stD.write(String.valueOf(stn) + ":" + String.valueOf(stRootedDis[stn]) + "\n");
        }
        stD.flush();
        stD.close();


        double totalD = 0.0;
        System.out.println("--------------------------------------------------------------");
        System.out.println("---------------------------Result-----------------------------");
        System.out.println("Running time:" + String.valueOf(costtime));
        System.out.println("Best Likelihood value : " + String.valueOf(GLOBALBESTLL));
        System.out.println("Best #iteration : " + bestLLIter);


        BufferedWriter ldOut = new BufferedWriter(new FileWriter(resultFolder + "distance.txt"));

        ldOut.write(astInit[0] +  ";" + "RF(AST_ST)" + "\n");
        ldOut.write(astInit[1]+ ";" + "RF(AST_GT)" + "\n");
        ldOut.write(astInit[2] +  ";" + "RF(GLASS_ST)" + "\n");
        ldOut.write(astInit[3] +  ";" + "RF(AST_ST by UPGMA)" + "\n");
        //ldOut.write(astInit[4] +  ";" + "RF(AST_WST by UPGMA)" + "\n");
        String stDistance = "!!!!!!RF of Best ST : " + operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
        ldOut.write(operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)) + ";RF(EM)" + "\n");
        System.out.println(stDistance);
        stDistance = "(Rooted_RF of Best ST) : " + operator.getRootDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
        ldOut.write(operator.getRootDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST)) + ";Rooted_RF(EM)" + "\n");
        System.out.println(stDistance);
        String gtDistance = "RF of GT_";
        double thisD = 0.0;
        for (int k = 0; k < lociNum; k++) {
            thisD = operator.getDistance(Trees.readTree(GLOBALBESTGTS.get(k)), Trees.readTree(trueGTS.get(k)));
            ldOut.write(gtDistance + k + " :" + thisD + "\n");
            totalD += thisD;
        }
        ldOut.write(totalD + "\n");
        ldOut.write(totalD / trueGTS.size() + "\n");
        System.out.println("Total distance of GTS is " + totalD + "\n");
        System.out.println("Average distance of GTS is " + totalD / trueGTS.size() + "\n");

        Network bestLocusST = inferSTByAST(locusBestGTS);
        System.out.println( "RF of Best RAxML trees' ST : " + operator.getDistance(Trees.readTree(bestLocusST.toString()), Trees.readTree(trueST)));
        ldOut.write("RF of Best RAxML trees' ST : " + operator.getDistance(Trees.readTree(bestLocusST.toString()), Trees.readTree(trueST))+ "\n");
        //ldOut.write("Best RAxML trees' ST : " + bestLocusST.toString()+ "\n");

        //bestLocusST = inferSTByWAST(locusBestGTS);
        //ldOut.write("RF of Best Weighted RAxML trees' ST : " + operator.getDistance(Trees.readTree(bestLocusST.toString()), Trees.readTree(trueST)));
        //ldOut.write("Rooted RF of Best Weighted RAxML trees' ST : " + operator.getRootDistance(Trees.readTree(bestLocusST.toString()), Trees.readTree(trueST))+ "\n");
        totalD = 0.0;
        for (int k = 0; k < lociNum; k++) {
            thisD = operator.getDistance(locusBestGTS.get(k), Trees.readTree(trueGTS.get(k)));
            ldOut.write(gtDistance + k + " :" + thisD + "\n");
            totalD += thisD;
        }
        ldOut.write("Total distance of RAxML best GTS is " + totalD + "\n");
        ldOut.write("Average distance of RAxML best GTS is " + totalD / trueGTS.size() + "\n");
        System.out.println("Total  distance of RAxML best GTS is " + totalD + "\n");
        System.out.println("Average distance of RAxML best GTS is " + totalD / trueGTS.size() + "\n");


        ldOut.write("---------------------------Global best ST change -----------------------------" + "\n");
        for (int i = 0; i < bestSTD.size(); i++) {
            ldOut.write(String.valueOf(bestITNum.get(i) + " : " + String.valueOf(bestSTD.get(i) + "  :  " + String.valueOf(bestSTLL.get(i)))) + "\n");
        }
        System.out.println("Average distance of GTS is " + totalD / trueGTS.size());
        System.out.println("Best ST is " + GLOBALBESTST);
//        llOut.write("Half AVE: " + operator.getAverage(gtsHLL) + "\n");
//        llOut.write("Half SD: " + operator.getStandardDevition(gtsHLL) + "\n");
//        llOut.write("Full AVE: " + operator.getAverage(gtsLL) + "\n");
//        llOut.write("Full SD: " + operator.getStandardDevition(gtsLL) + "\n");
//        llOut.write("Current AVE: " + operator.getAverage(currentITLL) + "\n");
//        llOut.write("Current SD: " + operator.getStandardDevition(currentITLL) + "\n");
        llOut.flush();
        llOut.close();
        ldOut.flush();
        ldOut.close();


        return finalGTS;
    }

    //  Input:  tempBestST
    //  Output: current best GTS
    public Network<NetNodeInfo> refineGeneSet(Network<NetNodeInfo> tempST) throws ParseException, IOException {
        List<Tree> currentTopo = simulator.simulateGeneTrees(tempST, null, refineGTSize); //operator.MS

        List<Tree> msTopos = currentTopo;
        //System.out.println("ms: " + msTopos.get(0));

        // operator.getDistinguishedTopo(currentTopo);
        if(GET_MS_DETAIL)
            operator.getMSDist(trueTopos, msTopos,msDist_Locus, itNum);

        List<Double> msTopoLL = operator.getGTSLLBySTYF(msTopos, tempST);
        System.out.println("msLL: " + msTopoLL.get(0));
        //TODO: involve p(gt|st)
        List<Tree> rGTS = getBestGTSByR(msTopos, tempST, msTopoLL);
        boolean move = false;
        double stDistance = 0.0;
        double gtDistance = 0.0;


        double p1 = currentHLL;
        double p2 = currentLL - currentHLL;
        currentITLL[itNum] = currentLL;
        currentITHLL[itNum] = currentHLL;
        currentITSTLL[itNum] = currentLL - currentHLL;
        double ifMove = (p1 - GLOBALBESTHLL) + ll_Ratio*(p2-(GLOBALBESTSTLL));
        if(itNum==1)
            ifMove = 1;
        if ( ifMove>0) {
            move = true;
            maxLL = currentLL;
            bestLLIter = itNum;
            GLOBALBESTST = tempST.toString();
            bestSTLL.add(currentLL);
            GLOBALBESTLL = currentLL;
            GLOBALBESTHLL = currentHLL;
            GLOBALBESTSTLL = GLOBALBESTLL-GLOBALBESTHLL;
            stDistance = operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
            bestSTD.add(stDistance);
            bestITNum.add(itNum);
            Iterator it = rGTS.iterator();
            GLOBALBESTGTS = new ArrayList<String>();
            while (it.hasNext()) {
                GLOBALBESTGTS.add(it.next().toString());
            }
            System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GLOBALBEST ST change");
            System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ST distance : " + stDistance);
        }
        if (STRATEGY == "NONE") {
            move = true;
        }
        else if (STRATEGY == "IMPROVE") {

        }
        else if (STRATEGY == "RANDOM") {
            if (currentLL >= maxLL) {
                move = true;
            }
            else {
                double random = Math.random();
                if (currentHLL > GLOBALBESTHLL) {
                    if (random < 0.95)
                        move = false;
                } else if ((p2 - GLOBALBESTSTLL) > 0) {
                    if (random < 0.95)
                        move = false;
                } else {
                    System.out.println("--------Random-Accept: #Iteration " + itNum + " : " + currentLL);
                    move = true;

                }
            }
        }
        if (move == true) {
            tempST = inferSTByAST(rGTS);
            maxLL = currentLL;
            stDis[itNum] = operator.getDistance(Trees.readTree(tempST.toString()), Trees.readTree(trueST));
            System.out.println("ST Distance of # " + itNum + " is : " + stDis[itNum]);
            stRootedDis[itNum] = operator.getRootDistance(Trees.readTree(tempST.toString()), Trees.readTree(trueST));
            System.out.println("Rooted ST Distance of # " + itNum + " is : " + stRootedDis[itNum]);
            double thisD = 0.0;
            for (int k = 0; k < lociNum; k++) {
                double temp = operator.getDistance(rGTS.get(k), Trees.readTree(trueGTS.get(k)));
                gtDis_Locus[itNum][k] = temp;
                thisD += temp;
            }
            gtDis[itNum] = thisD / lociNum;
            System.out.println("GTS Distance of # " + itNum + " is : " + gtDis[itNum]);
            System.out.println("**********Accept ST  ");
            gtsLL[itNum] = currentLL;
            gtsHLL[itNum] = currentHLL;
            endCounter = 0;

        } else {
            System.out.println("----------Reject MS GTS  ");
            gtsLL[itNum] = gtsLL[itNum - 1];
            gtsHLL[itNum] = gtsHLL[itNum - 1];

            stDis[itNum] = stDis[itNum - 1];
            stRootedDis[itNum] = stRootedDis[itNum - 1];
            System.out.println("ST Distance of # " + itNum + " is : " + stDis[itNum]);
            System.out.println("Rooted ST Distance of # " + itNum + " is : " + stRootedDis[itNum]);

            gtDis[itNum] = gtDis[itNum - 1];
            System.out.println("GTS Distance of # " + itNum + " is : " + gtDis[itNum]);

            for (int k = 0; k < lociNum; k++) {
                gtDis_Locus[itNum][k] = gtDis_Locus[itNum - 1][k];
            }
        }
        System.out.println("Current combined ll : " + currentLL);
        System.out.println("Full ll : " + gtsLL[itNum] + " ||| " + "Half ll : " + gtsHLL[itNum] + " ||| " + "Current best ll : " + GLOBALBESTLL);
        return tempST;
    }


    //Input: msTopos, ST, msTopoLL
    //Output: current GTS_nu no external, current LL
    public List<Tree> getBestGTSByR(List<Tree> topos, Network<NetNodeInfo> tempST, List<Double> msTopoLL) throws IOException {
        iterLLSeq = new double[lociNum];
        iterLLST = new double[lociNum];
        currentHLL = 0.0;
        currentLL = 0.0;
        List<Tree> bestGTs = new ArrayList<Tree>(); // scaled RAxML trees
        List<Double> bestGTLL = new ArrayList<Double>();
        List<Integer> bestGTNum = new ArrayList<Integer>();

        operator.updateMSTopos(topos);
        operator.runUpdateShell(refineGTSize, lociNum);
        for (int i = 0; i < lociNum; i++) {
            double[] llList = operator.getLocusLL(i, refineGTSize);

            double bestLL = llList[0] + msTopoLL.get(0);
            int bestGTNumi = 0;
            for (int j = 0; j < llList.length; j++) {
                if(llList[j]>locusGTLL[i]){
                    locusGTLL[i] = llList[j];
                    locusBestGTS.set(i,operator.getbestGTi(i, j));
                }

                double tempLL = llList[j] + msTopoLL.get(j);
                if (bestLL < tempLL) {
                    bestLL = tempLL;
                    bestGTNumi = j;
                }
            }
            bestGTs.add(operator.scaleGT(operator.getbestGTi(i, bestGTNumi), halfTheta, false));
            bestGTLL.add(bestLL);
            //TODO: get P(gt_j|st)
            currentLL +=  llList[bestGTNumi] + msTopoLL.get(bestGTNumi);
            bestGTNum.add(bestGTNumi);
            iterLLST[i] = msTopoLL.get(bestGTNumi);
            iterLLSeq[i] = llList[bestGTNumi];
            currentHLL += llList[bestGTNumi];
            //currentLL += bestLL;
        }

        return bestGTs;
    }


    // Use UPGMA + GLASS + AST
    // Input: Alignments
    // Output: AST ST + GLASS tree
    // GTS with BL: mutations/site
    public Network<NetNodeInfo> initST(List<Alignment> aln) throws IOException, ParseException, InterruptedException {

        //List<Tree> tempOGTs = operator.initByRAxML(aln, _seqLens[0],taxaNum+1, lociNum);
        List<UltrametricTree> uGTS = operator.simGTSByUPGMA(aln);

        operator.scaleGTS(uGTS, halfTheta, true);
        List<Tree> trees = operator.getTreeByUTree(uGTS);
        for(int i = 0; i<lociNum; i++){
            locusBestGTS.add(trees.get(i));
        }
        Constraint_ST = operator.inferSTByGLASS(trees);

        Network ast = inferSTByAST(trees);
        astInit[2] = operator.getDistance(Constraint_ST, Trees.readTree(trueST)); // (RF Glass)

        System.out.println("@@@@@AST init distance(unrooted): " + operator.getDistance(Trees.readTree(ast.toString()), Trees.readTree(trueST)));
        //   (RF AST)
        astInit[3] = operator.getDistance(Trees.readTree(ast.toString()), Trees.readTree(trueST));
        ast = inferSTByWAST(trees);

        INIT_ST = ast;
        // (RF wAST)
        astInit[4] =  operator.getDistance(Trees.readTree(ast.toString()), Trees.readTree(trueST));
        //TODO: add weight:
        //System.out.println("@@@@@WAST distance(unrooted): " + operator.getDistance(Trees.readTree(ast.toString()), Trees.readTree(trueST)));

        List<Double> gtsllList = operator.getGTSLLBySTYF(trees, Networks.readNetwork(Constraint_ST.toString()));
        currentHLL = Double.NEGATIVE_INFINITY;
        currentLL = 0.0;

//        for (int i = 0; i < gtsllList.size(); i++) {
//            currentLL += gtsllList.get(i);
//            System.out.println("currentGTLL : " + gtsllList.get(i));
//
//        }
        System.out.println("currentLL : " + currentLL);
        GLOBALBESTLL = Double.NEGATIVE_INFINITY;
        GLOBALBESTHLL = Double.NEGATIVE_INFINITY;
        GLOBALBESTSTLL =  Double.NEGATIVE_INFINITY;

        maxLL = Double.NEGATIVE_INFINITY;

        gtsHLL[0] = Double.NEGATIVE_INFINITY;
        gtsLL[0] = currentLL;
        GLOBALBESTST = INIT_ST.toString();
        // trueST = operator.rerootAndRemove(trueST,"O");
        stDis[0] = operator.getDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));
        stRootedDis[0] = operator.getRootDistance(Trees.readTree(GLOBALBESTST), Trees.readTree(trueST));

        double thisD = 0.0;
        //finalGTS = tempGTs;
        for (int k = 0; k < lociNum; k++) {
            //trueGTS.set(k,operator.rerootAndRemove(trueGTS.get(k),"O"));
            double d = operator.getDistance((Tree) trees.get(k), Trees.readTree(trueGTS.get(k)));
            System.out.println("# " + k + " : " + d);
            gtDis_Locus[0][k] = d;
            msDist_Locus[0][k] = 0;
            thisD += d;
        }
        gtDis[0] = thisD / lociNum;
        return INIT_ST;
    }


    //input: rooted gts with OG
    //output: rooted ST without OG
    public BniNetwork<NetNodeInfo> inferSTByAST(List<Tree> gts) throws IOException, ParseException {
        //String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i /Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";
        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i " + _ASTRALdir + "test_data/tempIn.tre -o " + _ASTRALdir + "test_data/testOut.tre"; ///Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";

        String cmdFile = _ASTRALdir + "tempCMD.sh";
        BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
        cmd.write(command + '\n');
        cmd.flush();
        cmd.close();
        String inputFile = _ASTRALdir + "test_data/tempIn.tre";
        BufferedWriter in = new BufferedWriter(new FileWriter(inputFile));
        Tree temp = gts.get(0);

        //TODO: add weight
        for (int i = 0; i < gts.size(); i++) {
            in.write(gts.get(i).toString() + "\n");
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
        for (int i = 1; i < stString.length(); i++) {
            if (stString.charAt(index1) == ':') {
                indexes.add(i);
            }
            index1++;
        }
        for (int i = 1; i <= indexes.size(); i++) {
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
        //System.out.println("AST String :" + stString);
        Tree st0 = Trees.readTree(stString);
        st0 = addExternal(st0, Constraint_ST);
        st0 = rootST(st0);
        BniNetwork n = (BniNetwork) Networks.readNetwork(st0.toString());
        return n;
    }

    //input: rooted gts with OG
    //output: rooted ST without OG
    public String initAST(List<Tree> gts) throws IOException, ParseException {
//        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i /Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";

        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i " + _ASTRALdir + "test_data/tempIn.tre -o " + _ASTRALdir + "test_data/testOut.tre"; ///Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";

        String cmdFile = _ASTRALdir + "tempCMD.sh";
        BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
        cmd.write(command + '\n');
        cmd.flush();
        cmd.close();
        String inputFile = _ASTRALdir + "test_data/tempIn.tre";
        BufferedWriter in = new BufferedWriter(new FileWriter(inputFile));
        Tree temp = gts.get(0);

        //TODO: add weight
        for (int i = 0; i < gts.size(); i++) {
            in.write(gts.get(i).toString() + "\n");
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
        for (int i = 1; i < stString.length(); i++) {
            if (stString.charAt(index1) == ':') {
                indexes.add(i);
            }
            index1++;
        }
        for (int i = 1; i <= indexes.size(); i++) {
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
        System.out.println("AST String :" + stString);

        return stString;
    }

    //input: rooted gts with OG
    //output: rooted ST without OG
    public BniNetwork<NetNodeInfo> inferSTByWAST(List<Tree> gts) throws IOException, ParseException {
        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i " + _ASTRALdir + "test_data/tempIn.tre -o " + _ASTRALdir + "test_data/testOut.tre"; ///Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";
        String cmdFile = _ASTRALdir + "tempCMD.sh";
        BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
        cmd.write(command + '\n');
        cmd.flush();
        cmd.close();
        String inputFile = _ASTRALdir + "test_data/tempIn.tre";
        BufferedWriter in = new BufferedWriter(new FileWriter(inputFile));
        Tree temp = gts.get(0);

        //TODO: add weight
        for (int i = 0; i < gts.size(); i++) {
            for(int j = (int)weights[i];j>0;j--){
                in.write(gts.get(i).toString() + "\n");

            }
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
        for (int i = 1; i < stString.length(); i++) {
            if (stString.charAt(index1) == ':') {
                indexes.add(i);
            }
            index1++;
        }
        for (int i = 1; i <= indexes.size(); i++) {
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
        System.out.println("AST String :" + stString);
        Tree st0 = Trees.readTree(stString);
        st0 = addExternal(st0, Constraint_ST);
        st0 = rootST(st0);
        BniNetwork n = (BniNetwork) Networks.readNetwork(st0.toString());
        return n;
    }


    //TODO
    //  Input: AST tree, glass tree
    //  Output: unrooted ST with external branches
    public Tree addExternal(Tree st0, Tree glassTree) throws IOException, ParseException {
        Network st = Networks.readNetwork(st0.toNewick());
        List<Tree> glassT = new ArrayList<Tree>();
        glassT.add(glassTree);
        Network newST = optimizeBL(glassT, st, 0.0001);
        st0 = Trees.readTree(newST.toString());
        return st0;
    }

    //  for all branches in st, if length = 0 -> length = ogHeight/10000
    //  if gts are ultrametric, add glassTree to bound ST
    public Network<NetNodeInfo> optimizeBL(List<Tree> gtsOG, Network st, Double ogHeight) {

        Map<Integer, List<String>> map = (Map<Integer, List<String>>) operator.getLevels(st);
        Set keys = map.keySet();
        int[] key = new int[keys.size()];
        Iterator keyIt = keys.iterator();
        for (int i = 0; i < keys.size(); i++) {
            key[i] = (int) keyIt.next();
        }
        int temp = 0;
        //bubble sort
        for (int i = 0; i < key.length - 1; i++) {
            for (int j = 0; j < key.length - 1; j++) {
                if (key[j] < key[j + 1]) {
                    temp = key[j];
                    key[j] = key[j + 1];
                    key[j + 1] = temp;
                }
            }
        }

        HashMap<String, Double> height = new HashMap<String, Double>();
        for (int i = 0; i < keys.size(); i++) {
            List<String> nodes = map.get(key[i]);

            Iterator<String> nodeIt = nodes.iterator();
            while (nodeIt.hasNext()) {
                BniNetNode tempNode = (BniNetNode) st.findNode(nodeIt.next());
                BniNetNode p = (BniNetNode) tempNode.getParents().iterator().next();
                if (tempNode.getParentDistance(p) != Double.NEGATIVE_INFINITY)
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
                        if (distance == 0 || distance == Double.POSITIVE_INFINITY) {
                            distance = ogHeight;
                        }
                        tempNode.setParentDistance(p, distance / 2);
                        sib.setParentDistance(p, distance / 2);
                        height.put(p.getName(), distance / 2);

                    }
                    //sib is non-leaf
                    else {
                        // lowerDist
                        if (sib.getParentDistance(p) == 0 || sib.getParentDistance(p) == Double.NEGATIVE_INFINITY)
                            //TODO
                            sib.setParentDistance(p, ogHeight / 100);
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
        if (distance == Double.POSITIVE_INFINITY)
            System.err.println("Wrong Distance!!" + nodeName1 + " + " + nodeName2);

        return distance;
    }


    //TODO
    //  Input: unrooted AST tree with external branches
    //  Output: rooted & ultrametric ST
    public Tree rootST(Tree st0) throws IOException, ParseException {
        //get tree height
        STITree rootT = (STITree) st0;
        Trees.autoLabelNodes(rootT);
        TNode root = rootT.getRoot();
        Iterator subT = root.getChildren().iterator();
        TNode left = (TNode) subT.next();
        TNode right = (TNode) subT.next();
        String[] leftStr = getNodeHeight(left);
        String[] rightStr = getNodeHeight(right);
        double h1 = Double.valueOf(leftStr[1]);
        double h2 = Double.valueOf(rightStr[1]);
        Double treeHeight = ( h1 + h2 )/2;
        String start = (h1 > h2)? leftStr[0]: rightStr[0];
        TNode sNode = rootT.getNode(start);

        //find root edge
        double dist = 0.0;
        STINode newRoot = (STINode) sNode;// new root
        TNode fRoot = sNode;// root edge finish node
        while(true){
            dist += sNode.getParentDistance();
            if(dist >= treeHeight){
                //no change node
                if(dist == treeHeight){
                    newRoot = (STINode)sNode.getParent();
                    rootT.rerootTreeAtNode(newRoot);
                }
                else if(dist > treeHeight) {
                    double diff = dist - treeHeight;
                    double edgeLength = sNode.getParentDistance();
                    BniNetwork net = (BniNetwork) Networks.readNetwork(rootT.toString());
                    BniNetNode newNetRoot = new BniNetNode("I" + String.valueOf(rootT.getNodeCount()),0.0);
                    BniNetNode oldP = (BniNetNode)net.findNode(sNode.getParent().getName());
                    oldP.adoptChild(newNetRoot,diff);
                    newNetRoot.adoptChild((BniNetNode)net.findNode(sNode.getName()),edgeLength-diff);
                    oldP.removeChild(net.findNode(sNode.getName()));

                    BniNetNode oldNetRoot = (BniNetNode)net.getRoot();
                    BniNetNode tempC = newNetRoot;
                    BniNetNode tempP = oldP;
                    while(!tempC.equals(oldNetRoot)){
                        double tempD = tempC.getParentDistance(tempP);
                        tempP.removeChild(tempC);
                        tempC.adoptChild(tempP,tempD);
                        tempC = tempP;
                        tempP = (BniNetNode) tempC.getParents().iterator().next();
                    }
                    BniNetNode fChild = (BniNetNode) tempC.getChildren().iterator().next();
                    tempP.removeChild(tempC);
                    tempP.adoptChild(fChild,fChild.getParentDistance(tempC));
                    tempC.removeChild(fChild);
                    net.resetRoot(newNetRoot);
                    rootT = (STITree) Trees.readTree(net.toString());
                }
                break;
            }
            else{
                sNode = sNode.getParent();
            }
        }
        Tree result = operator.getUTree(rootT,0.0);
        return result;
    }
//TODO keep best gts P(seq|GT)


    //Input: Node name
    //Output: String[2]   0:leaf name  1： leaf distance
    public String[] getNodeHeight(TNode temp) {
        if (temp.isLeaf()){
            String[] str = new String[2];
            str[0] = temp.getName();
            str[1] = String.valueOf(0.0);
            return str;
        }
        else {
            String[] tempStr1;
            String[] tempStr2;
            String[] tempResult = new String[2];
            Iterator it = temp.getChildren().iterator();
            TNode left = (TNode) it.next();
            tempStr1 = getNodeHeight(left);
            double leftD = Double.valueOf(tempStr1[1]) + left.getParentDistance();
            TNode right = (TNode) it.next();
            tempStr2 = getNodeHeight(right);
            double rightD = Double.valueOf(tempStr2[1]) + right.getParentDistance();
            if(leftD>rightD){
                tempResult[0] = tempStr1[0];
                tempResult[1] = String.valueOf(leftD);
            }
            else{
                tempResult[0] = tempStr2[0];
                tempResult[1] = String.valueOf(rightD);
            }
            return tempResult;
        }
    }



}

