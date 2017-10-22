package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import sun.nio.ch.Net;

import java.io.*;
import java.util.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by doriswang on 3/12/17.
 */

public class IIGTSimulator {
    //TODO: change path
    protected static int numGTs;
    //private static double[] _scales = new double[] {0.1, 0.25, 0.5, 1.0};
    private static double[] _scales;// = new double[] {0.5};
    private static int[] _seqLens;// = new int[] {500};
    //private static int[] _seqLens = new int[] {250, 500, 1000};
    //private static int[] _numLocis = new int[] {32, 64, 128};
    private static double halfTheta;
    private static int iteration;

    private static String msPath = "/Users/doriswang/PhyloNet/tools/msFiles/msdir/ms";
    private static String seqGenePath = "/Users/doriswang/PhyloNet/tools/Seq-Gen.v1.3.3/source/seq-gen";
    private static String TREE_DIR = "/Users/doriswang/PhyloNet/Data/IIG/tree/";
    private static String SEQ_DIR = "/Users/doriswang/PhyloNet/Data/IIG/seq/";
    private static String _DIR =
            "/Users/doriswang/PhyloNet/Data/IIG/";

    IIGTSimulator(int numGT, double[] scales, int[] seqLens, double halfTheta, int t) {
        this.numGTs = numGT;
        this._scales = scales;
        this._seqLens = seqLens;
        this.halfTheta = halfTheta;
        this.iteration = t;
    }

    public static void main(String[] args) throws IOException, InterruptedException{
        _seqLens = new int[] {200};
        IIGTSimulator simulator1 = new IIGTSimulator(100, _scales, _seqLens, 0.0005, 10);
        simulator1.simulateSampleData();
//        simulator1.simulateSeqByGTS(32, 1000, 0.005, "/Users/doriswang/PhyloNet/Data/17-taxon/001/ST0/1/Seq/");

            //RawDataProcessor r = new RawDataProcessor();
        //r.compare();

//        Tree t = Trees.readTree("(1:0.03798,(((((2:0.004797,9:0.003204,3:0.000806,4:0.000772,5:0.000810,6:0.000787,7:0.000786,10:0.000798):0.001607,8:0.002159):0.037802,(((11:0.000817,13:0.000808,14:0.000809):0.003299,(12:0.000806,15:0.000804,16:0.000819,17:0.000806):0.002343):0.002431,(18:0.000788,19:0.000807):0.001933):0.029201):0.032477,((20:0.000816,(21:0.000801,22:0.000811,23:0.000801):0.004790,24:0.000806):0.001593,25:0.001260):0.049365):0.076642):0.031080);");
//        TNode n = t.getRoot();
//        MutableTree mT = (MutableTree)t;
//        Trees.autoLabelNodes(mT);
//        System.out.println(mT.toString());
//        InferOperator op = new InferOperator(10);
//        double d = op.getDistance(t,(Tree)mT);
//        System.out.println(d);

    }

    //Input: gts, seqLen, halfTheta, seqPath
    //Output: seq files
    public static List<Tree> simulateSeqByGTS(int lociNum, int seqLen, double htheta, String seqPath) throws IOException, InterruptedException{
        List<Tree> gTrees = new ArrayList<Tree>();
        Network trueST = Networks.readNetwork("(O:4000000.0000000005,((7:320410.56025586306,(5:171859.33416363757,2:171859.33416363757):148551.22609222552):479589.4397441369,((((15:1286.8750574758267,1:1286.8750574758267):193151.5628134453,(6:190005.17240669002,(10:31789.288803547137,3:31789.288803547137):158215.88360314287):4433.265464231158):70363.8723491511,(((8:47699.30894420634,12:47699.30894420634):19910.799516858293,11:67610.10846106464):18959.454862388262,(4:44421.322659187055,16:44421.322659187055):42148.240664265846):178232.74689661936):478658.742741006,((9:104854.84644501805,14:104854.84644501805):36171.64602825128,13:141026.4924732693):602434.560487809):56538.947038921746):3200000.0000000005);");
        Iterator itNode = trueST.getTreeNodes().iterator();
        while (itNode.hasNext()) {
            NetNode n = (NetNode) itNode.next();
            if (n == (NetNode) trueST.getRoot())
                continue;

            NetNode p = (NetNode) n.getParents().iterator().next();
            //double ratio = currentHeight - heightN;
            if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                n.setParentDistance(p, n.getParentDistance(p) / 1000000); // scale to c_Unit
            }
        }
        InferOperator ifo = new InferOperator(1000);
        ifo.isExitsPath(seqPath);
        BufferedWriter stbw = new BufferedWriter(new FileWriter(seqPath + "trueST.txt"));
        stbw.write(trueST.toString());
        stbw.flush();
        stbw.close();

        BufferedReader gtReader = new BufferedReader(new FileReader("/Users/doriswang/PhyloNet/Data/17-taxon/32loci/Rep0gtrees"));

        //String folder = "" + iteration;
        String seqDir = "/Users/doriswang/PhyloNet/Data/17-taxon/001/ST0/1/Seq/";
        String treeDir = seqPath + "gtree/";
        ifo.isExitsPath(treeDir);

        BufferedWriter bw = new BufferedWriter(new FileWriter(seqPath + "trueGTS.txt"));

        for (int tid = 0; tid < lociNum; tid++) {

            String gt = gtReader.readLine().trim();
            Tree tempGT = Trees.readTree(gt);
            Iterator gtNode = tempGT.getNodes().iterator();
            while (gtNode.hasNext()) {
                TNode n = (TNode) gtNode.next();
                if (n == (TNode) tempGT.getRoot())
                    continue;

                TNode p = (TNode) n.getParent();
                //double ratio = currentHeight - heightN;
                if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(n.getParentDistance() / 1000000); // scale to c_Unit
                }
            }
            bw.write(tempGT.toString() + "\n");
            try {
                String fileName = treeDir + "TrueGT_" + tid + ".tree";
                BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                out.write(tempGT.toString() + "\n");
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            String command = seqGenePath + " -mHKY -s" + htheta + " -l "
                    + seqLen + " -on < " + treeDir + "TrueGT_" + tid + ".tree > "
                    + seqPath + "seq_" + tid + "_" + seqLen + ".nex";
            //SEQ_DIR + "/" + iteration + "/seq_" + seqNum + "_" + seqLength +".nex";
            try {
                String fileName = seqPath + "seqgen_" + tid + "_" + seqLen + ".sh";
                BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                out.write(command + "\n");
                out.close();
                ProcessBuilder pb = new ProcessBuilder("/bin/bash", fileName);
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
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        bw.close();
        return gTrees;
    }


    public static List<String> tempSimulateData() throws IOException{
        List<String> trees = new ArrayList<String>();
                Network trueST = Networks.readNetwork("(O:4000000.0000000005,((7:320410.56025586306,(5:171859.33416363757,2:171859.33416363757):148551.22609222552):479589.4397441369,((((15:1286.8750574758267,1:1286.8750574758267):193151.5628134453,(6:190005.17240669002,(10:31789.288803547137,3:31789.288803547137):158215.88360314287):4433.265464231158):70363.8723491511,(((8:47699.30894420634,12:47699.30894420634):19910.799516858293,11:67610.10846106464):18959.454862388262,(4:44421.322659187055,16:44421.322659187055):42148.240664265846):178232.74689661936):478658.742741006,((9:104854.84644501805,14:104854.84644501805):36171.64602825128,13:141026.4924732693):602434.560487809):56538.947038921746):3200000.0000000005);");
                Iterator itNode = trueST.getTreeNodes().iterator();
                while (itNode.hasNext()) {
                    NetNode n = (NetNode) itNode.next();
                    if (n == (NetNode) trueST.getRoot())
                        continue;

                    NetNode p = (NetNode) n.getParents().iterator().next();
                    //double ratio = currentHeight - heightN;
                    if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(p, n.getParentDistance(p) / 40000); // Coalescent_Unit = Tao /  theta/2
                    }
                }
                trees.add(trueST.toString());
                //String folder = "" + iteration;
                String seqDir = "/Users/doriswang/PhyloNet/Data/17-taxon/004/ST0/3/Seq/";
                String treeDir = "/Users/doriswang/PhyloNet/Data/17-taxon/004/ST0/3/Tree/";

                BufferedWriter bw = new BufferedWriter(new FileWriter(treeDir + "trueST.txt"));
                bw.write(trueST.toString() + "\n");
                Map<String, List<String>> species2alleles = null;
                List<Tree> trueGTs = simulateGeneTrees(trueST, species2alleles, 10);
                for (int tid = 0; tid < trueGTs.size(); tid++) {
                    String gt = trueGTs.get(tid).toNewick();
                    trees.add(gt);
                    bw.write(gt + "\n");
                    try {
                        String fileName = treeDir + "TrueGT_" + tid + ".tree";
                        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                        out.write(gt + "\n");
                        out.close();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    for (int seqLen : _seqLens) {
                        String command = seqGenePath + " -mHKY -s" + halfTheta + " -l "
                                + seqLen + " -on < " + treeDir + "TrueGT_" + tid + ".tree > "
                                + seqDir + "seq_" + tid + "_" + seqLen + ".nex";
                        //SEQ_DIR + "/" + iteration + "/seq_" + seqNum + "_" + seqLength +".nex";
                        try {
                            String fileName = seqDir + "seqgen_" + tid + "_" + seqLen + ".sh";
                            BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                            out.write(command + "\n");
                            out.close();
                            ProcessBuilder pb = new ProcessBuilder("/bin/bash", fileName);
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
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }
                bw.close();

        return trees;
    }

    public static List<String> loadTrees( String filePath, int lociNum) throws IOException{
        List<String> trees = new ArrayList<String>();
        String streeFile = filePath + "Tree/trueST.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
        String trueST = (String) stReader.readLine().trim();
        trees.add(trueST);
        for (int ln = 0; ln < lociNum ; ln++) {
            String tree = stReader.readLine().trim();
            trees.add(tree);
        }
        return trees;
    }


    public static List<String> simulateData() {
        List<String> trees = new ArrayList<String>();
        try {
            //para: theta(mutation_rate)  = 4 * pop_size * mu
            //double halfTheta = 0.008;
            //THETA = 0.016

            for (double scale : _scales) {


                //Network net =  getNetWithO(scale);
                //para: trueST-> true species tree with branch
                Network trueST = get16ST(scale);
                System.out.println(trueST.toString());

                //For 17 Taxon
//                String st = "((C:0.002859330148428049,G:0.002859330148428049):0.010473421735219215,(R:0.010936563366414687,(L:0.008340943025951807,(A:0.004491653225077681,Q:0.004491653225077681):0.0038492898008741254):0.00259562034046288):0.0023961885172325784);";
//                Network trueST =  Networks.readNetwork(st);
                Iterator itNode = trueST.getTreeNodes().iterator();
                while (itNode.hasNext()) {
                    NetNode n = (NetNode) itNode.next();
                    if (n == (NetNode) trueST.getRoot())
                        continue;

                    NetNode p = (NetNode) n.getParents().iterator().next();
                    //double ratio = currentHeight - heightN;
                    if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(p, n.getParentDistance(p) / halfTheta); // Coalescent_Unit = Tao /  theta/2
                    }
                }

                trees.add(trueST.toString());


                //String folder = "" + iteration;
                String seqDir = SEQ_DIR + iteration + "/";
                String treeDir = TREE_DIR + iteration + "/";

                BufferedWriter bw = new BufferedWriter(new FileWriter(treeDir + "trueST_" + iteration + ".txt"));
                bw.write(trueST.toString() + "\n");
                Map<String, List<String>> species2alleles = null;
                List<Tree> trueGTs = simulateGeneTrees(trueST, species2alleles, numGTs);
                for (int tid = 0; tid < trueGTs.size(); tid++) {
                    String gt = trueGTs.get(tid).toNewick();
                    trees.add(gt);
                    bw.write(gt + "\n");
                    try {
                        String fileName = treeDir + "TrueGT_" + tid + ".tree";
                        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                        out.write(gt + "\n");
                        out.close();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    for (int seqLen : _seqLens) {
                        String command = seqGenePath + " -mHKY -s" + halfTheta + " -l "
                                + seqLen + " -on < " + treeDir + "TrueGT_" + tid + ".tree > "
                                + seqDir + "seq_" + tid + "_" + seqLen + ".nex";
                        //SEQ_DIR + "/" + iteration + "/seq_" + seqNum + "_" + seqLength +".nex";
                        try {
                            String fileName = seqDir + "seqgen_" + tid + "_" + seqLen + ".sh";
                            BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                            out.write(command + "\n");
                            out.close();
                            ProcessBuilder pb = new ProcessBuilder("/bin/bash", fileName);
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
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }
                bw.close();

//                Network outputNet = getTrueST(scale);
//                Networks.scaleNetwork(outputNet, halfTheta);
//
//                BufferedWriter bw2 = new BufferedWriter(new FileWriter(seqDir + "startNetGTs_" + iteration + ".txt"));
//                bw2.write("[" + (halfTheta * 2) + "]" + outputNet.toString() + "\n");
//
//                for(Tree t : trueGTs) {
//                    //Trees.rootAndRemoveOutgroup((STITree) t, "O");
//                    Trees.scaleBranchLengths((STITree) t, halfTheta);
//                    bw2.write(t.toString() + "\n");
//                }
//                bw2.close();

            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return trees;
    }

    public static List<String> simulateSampleData() {
        List<String> trees = new ArrayList<String>();
        try {
                Network trueST = Networks.readNetwork("(((((1:1,2:1):1,(3:1,4:1):1):1,((5:1,6:1):1,(7:1,8:1):1):1):1,(((9:1,10:1):1,(11:1,12:1):1):1,((13:1,14:1):1,(15:1,16:1):1):1):1):16,O:20);");

                Iterator itNode = trueST.getTreeNodes().iterator();
                while (itNode.hasNext()) {
                    NetNode n = (NetNode) itNode.next();
                    if (n == (NetNode) trueST.getRoot())
                        continue;

                    NetNode p = (NetNode) n.getParents().iterator().next();
                    //double ratio = currentHeight - heightN;
                    if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(p, n.getParentDistance(p) ); // Coalescent_Unit = Tao /  theta/2
                    }
                }

                trees.add(trueST.toString());


                //String folder = "" + iteration;
                String seqDir = "/Users/doriswang/PhyloNet/Data/IIG/symmetrical/" + _seqLens[0] + "/";
                String treeDir = "/Users/doriswang/PhyloNet/Data/IIG/symmetrical/"  + _seqLens[0] + "/Tree/";

                BufferedWriter bw = new BufferedWriter(new FileWriter(treeDir + "trueST.txt"));
                bw.write(trueST.toString() + "\n");
                Map<String, List<String>> species2alleles = null;
                List<Tree> trueGTs = simulateGeneTrees(trueST, species2alleles, 100);
                for (int tid = 0; tid < 100; tid++) {
                    String gt = trueGTs.get(tid).toNewick();
                    trees.add(gt);
                    bw.write(gt + "\n");
                    try {
                        String fileName = treeDir + "TrueGT_" + tid + ".tree";
                        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                        out.write(gt + "\n");
                        out.close();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    for (int seqLen : _seqLens) {
                        String command = seqGenePath + " -mHKY -s" + halfTheta + " -l "
                                + seqLen + " -on < " + treeDir + "TrueGT_" + tid + ".tree > "
                                + seqDir + "seq_" + tid + "_" + seqLen + ".nex";
                        //SEQ_DIR + "/" + iteration + "/seq_" + seqNum + "_" + seqLength +".nex";
                        try {
                            String fileName = seqDir + "seqgen_" + tid + "_" + seqLen + ".sh";
                            BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
                            out.write(command + "\n");
                            out.close();
                            ProcessBuilder pb = new ProcessBuilder("/bin/bash", fileName);
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
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }
                bw.close();

//                Network outputNet = getTrueST(scale);
//                Networks.scaleNetwork(outputNet, halfTheta);
//
//                BufferedWriter bw2 = new BufferedWriter(new FileWriter(seqDir + "startNetGTs_" + iteration + ".txt"));
//                bw2.write("[" + (halfTheta * 2) + "]" + outputNet.toString() + "\n");
//
//                for(Tree t : trueGTs) {
//                    //Trees.rootAndRemoveOutgroup((STITree) t, "O");
//                    Trees.scaleBranchLengths((STITree) t, halfTheta);
//                    bw2.write(t.toString() + "\n");
//                }
//                bw2.close();


        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return trees;
    }


    //INITIAL SPECIES TREE SETTING
    public static Network getTrueST(double scale) {
        Network<Double> net = new BniNetwork<>();
        net.createRoot("I1");
        NetNode<Double> I1 = net.getRoot();
        NetNode<Double> I2 = new BniNetNode<Double>("I2", null);
        NetNode<Double> I3 = new BniNetNode<Double>("I3", null);
        NetNode<Double> I4 = new BniNetNode<Double>("I4", null);
        //CHANGE from         NetNode<Double> A = new BniNetNode<>("A", null);
        NetNode<Double> A = new BniNetNode<Double>("A", null);
        NetNode<Double> C = new BniNetNode<Double>("C", null);
        NetNode<Double> Q = new BniNetNode<Double>("Q", null);
        NetNode<Double> R = new BniNetNode<Double>("R", null);
        NetNode<Double> L = new BniNetNode<Double>("L", null);

        I1.adoptChild(I2, 0.5 * scale);
        I1.adoptChild(I3, 1.0 * scale);
        //I5.adoptChild(R, 1.0 * scale);
        //I5.adoptChild(I6, 0.25 * scale);
        I2.adoptChild(L, 1.5 * scale);
        I2.adoptChild(I4, 1.0 * scale);
        I3.adoptChild(Q, 1.0 * scale);
        I3.adoptChild(R, 1.0 * scale);
        //I6.adoptChild(Q, 0.75 * scale);
        I4.adoptChild(A, 0.5 * scale);
        I4.adoptChild(C, 0.5 * scale);

        //I6.setParentProbability(I5, 0.3);
        //I6.setParentProbability(I3, 0.7);

        return net;
    }

    //TODO:6
//    public static Network get16ST(double scale) {
//        Networks.readNetwork(())
//        return net;
//    }

    public static Network get16ST(double scale) {
        Network<Double> net = new BniNetwork<>();
        net.createRoot("I0");
        NetNode<Double> I0 = net.getRoot();
        NetNode<Double> I1 = new BniNetNode<Double>("I1", null);
        NetNode<Double> I2 = new BniNetNode<Double>("I2", null);
        NetNode<Double> I3 = new BniNetNode<Double>("I3", null);
        NetNode<Double> I4 = new BniNetNode<Double>("I4", null);
        NetNode<Double> I5 = new BniNetNode<Double>("I5", null);
        NetNode<Double> I6 = new BniNetNode<Double>("I6", null);
        NetNode<Double> I7 = new BniNetNode<Double>("I7", null);
        NetNode<Double> I8 = new BniNetNode<Double>("I8", null);
        NetNode<Double> I9 = new BniNetNode<Double>("I9", null);
        NetNode<Double> I10 = new BniNetNode<Double>("I10", null);
        NetNode<Double> I11 = new BniNetNode<Double>("I11", null);
        NetNode<Double> I12 = new BniNetNode<Double>("I12", null);
        NetNode<Double> I13 = new BniNetNode<Double>("I13", null);
        NetNode<Double> I14 = new BniNetNode<Double>("I14", null);
        //NetNode<Double> I15 = new BniNetNode<Double>("I7", null);
        //NetNode<Double> I8 = new BniNetNode<Double>("I8", null);

        NetNode<Double> A = new BniNetNode<Double>("A", null);
        NetNode<Double> B = new BniNetNode<Double>("B", null);
        NetNode<Double> C = new BniNetNode<Double>("C", null);
        NetNode<Double> D = new BniNetNode<Double>("D", null);
        NetNode<Double> E = new BniNetNode<Double>("E", null);
        NetNode<Double> F = new BniNetNode<Double>("F", null);
        NetNode<Double> G = new BniNetNode<Double>("G", null);
        NetNode<Double> H = new BniNetNode<Double>("H", null);
        NetNode<Double> I = new BniNetNode<Double>("I", null);
        NetNode<Double> J = new BniNetNode<Double>("J", null);
        NetNode<Double> K = new BniNetNode<Double>("K", null);
        NetNode<Double> L = new BniNetNode<Double>("L", null);
        NetNode<Double> M = new BniNetNode<Double>("M", null);
        NetNode<Double> N = new BniNetNode<Double>("N", null);
        NetNode<Double> O = new BniNetNode<Double>("O", null);
        NetNode<Double> P = new BniNetNode<Double>("P", null);

        I0.adoptChild(I1, 3.0);
        I0.adoptChild(I2, 1.0);

        I1.adoptChild(I3, 1.0);
        I1.adoptChild(P, 3.0);

        I2.adoptChild(I5, 4.0);
        I2.adoptChild(I4, 1.0);

        I3.adoptChild(N, 2.0);
        I3.adoptChild(O, 2.0);

        I4.adoptChild(I7, 1.0);
        I4.adoptChild(I6, 1.0);

        I5.adoptChild(A, 1.0);
        I5.adoptChild(I8, 0.5);


        I6.adoptChild(I10, 1.0);
        I6.adoptChild(I9, 2.0);

        I7.adoptChild(I12, 2.0);
        I7.adoptChild(I11, 1.0);

        I8.adoptChild(B, 0.5);
        I8.adoptChild(C, 0.5);

        I9.adoptChild(L, 1.0);
        I9.adoptChild(M, 1.0);

        I10.adoptChild(I13, 1.5);
        I10.adoptChild(K, 2.0);

        I11.adoptChild(F, 2.0);
        I11.adoptChild(I14, 1.0);

        I12.adoptChild(D, 1.0);
        I12.adoptChild(E, 1.0);

        I13.adoptChild(I, 0.5);
        I13.adoptChild(J, 0.5);

        I14.adoptChild(G, 1.0);
        I14.adoptChild(H, 1.0);

        return net;
    }


    //10 nodes: (((J,I)I6,(H,(G,F)I8)I5)I2,((E,D)I4,((C,B)I7,A)I3)I1)I0;
    //INITIAL SPECIES TREE SETTING
    public static Network getTrueBigST(double scale) {
        Network<Double> net = new BniNetwork<>();
        net.createRoot("I0");
        NetNode<Double> I0 = net.getRoot();
        NetNode<Double> I1 = new BniNetNode<Double>("I1", null);
        NetNode<Double> I2 = new BniNetNode<Double>("I2", null);
        NetNode<Double> I3 = new BniNetNode<Double>("I3", null);
        NetNode<Double> I4 = new BniNetNode<Double>("I4", null);
        NetNode<Double> I5 = new BniNetNode<Double>("I5", null);
        NetNode<Double> I6 = new BniNetNode<Double>("I6", null);
        NetNode<Double> I7 = new BniNetNode<Double>("I7", null);
        NetNode<Double> I8 = new BniNetNode<Double>("I8", null);

        NetNode<Double> A = new BniNetNode<Double>("A", null);
        NetNode<Double> B = new BniNetNode<Double>("B", null);
        NetNode<Double> C = new BniNetNode<Double>("C", null);
        NetNode<Double> D = new BniNetNode<Double>("D", null);
        NetNode<Double> E = new BniNetNode<Double>("E", null);
        NetNode<Double> F = new BniNetNode<Double>("F", null);
        NetNode<Double> G = new BniNetNode<Double>("G", null);
        NetNode<Double> H = new BniNetNode<Double>("H", null);
        NetNode<Double> I = new BniNetNode<Double>("I", null);
        NetNode<Double> J = new BniNetNode<Double>("J", null);

        I0.adoptChild(I1, 2.0);
        I0.adoptChild(I2, 0.5);

        I1.adoptChild(A, 1.0);
        I1.adoptChild(I3, 0.5);

        I2.adoptChild(I4, 1.0);
        I2.adoptChild(I5, 0.5);

        I3.adoptChild(B, 0.5);
        I3.adoptChild(C, 0.5);

        I4.adoptChild(I6, 0.5);
        I4.adoptChild(F, 1.5);

        I5.adoptChild(I7, 0.5);
        I5.adoptChild(I8, 1.0);

        I6.adoptChild(D, 1.0);
        I6.adoptChild(E, 1.0);

        I7.adoptChild(G, 1.5);
        I7.adoptChild(H, 1.5);

        I8.adoptChild(I, 1.0);
        I8.adoptChild(J, 1.0);


        return net;
    }

    //10 nodes: (((J,I)I6,(H,(G,F)I8)I5)I2,((E,D)I4,((C,B)I7,A)I3)I1)I0;
    //INITIAL SPECIES TREE SETTING
    public static Network get10ST(double scale) {
        Network<Double> net = new BniNetwork<>();
        net.createRoot("I1");
        NetNode<Double> I0 = net.getRoot();
        NetNode<Double> I1 = new BniNetNode<Double>("I1", null);
        NetNode<Double> I2 = new BniNetNode<Double>("I2", null);
        NetNode<Double> I3 = new BniNetNode<Double>("I3", null);
        NetNode<Double> I4 = new BniNetNode<Double>("I4", null);
        NetNode<Double> I5 = new BniNetNode<Double>("I5", null);
        NetNode<Double> I6 = new BniNetNode<Double>("I6", null);
        NetNode<Double> I7 = new BniNetNode<Double>("I7", null);
        NetNode<Double> I8 = new BniNetNode<Double>("I8", null);

        NetNode<Double> A = new BniNetNode<Double>("A", null);
        NetNode<Double> B = new BniNetNode<Double>("B", null);
        NetNode<Double> C = new BniNetNode<Double>("C", null);
        NetNode<Double> D = new BniNetNode<Double>("D", null);
        NetNode<Double> E = new BniNetNode<Double>("E", null);
        NetNode<Double> F = new BniNetNode<Double>("F", null);
        NetNode<Double> G = new BniNetNode<Double>("G", null);
        NetNode<Double> H = new BniNetNode<Double>("H", null);
        NetNode<Double> I = new BniNetNode<Double>("I", null);
        NetNode<Double> J = new BniNetNode<Double>("J", null);

        I7.adoptChild(B, 0.5);
        I7.adoptChild(C, 0.5);

        I3.adoptChild(A, 1.0);
        I3.adoptChild(I7, 0.5);

        I4.adoptChild(D, 1.0);
        I4.adoptChild(E, 1.0);

        I1.adoptChild(I3, 1.0);
        I1.adoptChild(I4, 1.0);

        I8.adoptChild(F, 0.5);
        I8.adoptChild(G, 0.5);

        I5.adoptChild(I8, 0.5);
        I5.adoptChild(H, 1.0);

        I6.adoptChild(I, 1.0);
        I6.adoptChild(J, 1.0);

        I2.adoptChild(I5, 1.0);
        I2.adoptChild(I6, 1.0);

        I0.adoptChild(I1, 1.0);
        I0.adoptChild(I2, 1.0);

        return net;
    }

    //INITIAL SPECIES TREE SETTING WITH OUTGROUP
    private static Network getTrueSTWithO(double scale) {
        Network<Double> net = new BniNetwork<>();
        net.createRoot("I0");
        NetNode<Double> I0 = net.getRoot();
        NetNode<Double> I1 = new BniNetNode<Double>("I1", null);
        NetNode<Double> I2 = new BniNetNode<Double>("I2", null);
        NetNode<Double> I3 = new BniNetNode<Double>("I3", null);
        NetNode<Double> I4 = new BniNetNode<Double>("I4", null);
        //CHANGE from         NetNode<Double> A = new BniNetNode<>("A", null);
        NetNode<Double> A = new BniNetNode<Double>("A", null);
        NetNode<Double> C = new BniNetNode<Double>("C", null);
        NetNode<Double> Q = new BniNetNode<Double>("Q", null);
        NetNode<Double> R = new BniNetNode<Double>("R", null);
        NetNode<Double> L = new BniNetNode<Double>("L", null);
        NetNode<Double> O = new BniNetNode<Double>("O", null);

        I0.adoptChild(I1, 10.0 - 2 * scale);
        I0.adoptChild(O, 10.0);
        I1.adoptChild(I2, 0.5 * scale);
        I1.adoptChild(I3, 1.0 * scale);
        //I5.adoptChild(R, 1.0 * scale);
        //I5.adoptChild(I6, 0.25 * scale);
        I2.adoptChild(L, 1.5 * scale);
        I2.adoptChild(I4, 1.0 * scale);
        I3.adoptChild(Q, 1.0 * scale);
        I3.adoptChild(R, 1.0 * scale);
        //I6.adoptChild(Q, 0.75 * scale);
        I4.adoptChild(A, 0.5 * scale);
        I4.adoptChild(C, 0.5 * scale);

        //I6.setParentProbability(I5, 0.3);
        //I6.setParentProbability(I3, 0.7);

        return net;
    }


    public static List<Tree> simulateGeneTrees(Network net,
                                               Map<String, List<String>> species2alleles, int numGT) {
        SimGTInNUNetByMS sim = new SimGTInNUNetByMS();
        List<Tree> gts = sim.generateGTs(Networks.readNetwork(net.toString()), species2alleles, numGT, msPath);
        //System.out.println(sim.getMSCommand());
        return gts;
    }

//    //nodesName can be used if leafName != int number
//    public static String generateMSCMD(Network net, int numGT, Map<String, Integer> nodesName) {
//        //List<Tree> gts =  new ArrayList<Tree>();
//        //net without outGroup
//        Networks.autoLabelNodes(net);
//        String msCMD = " " + net.getLeafCount() + " " + numGT + " -T -I" + " " + net.getLeafCount();
//        for (int i = 0; i < net.getLeafCount(); i++) {
//            msCMD += " 1";
//        }
//        System.out.println("MSCMD : " + msCMD.toString());
//
//        Map dist2Root = new HashMap<String, Double>();
//        Map heights = new HashMap<String, Double>();
//        Map in2leaf = new HashMap<String,String>();  //name mapping
//        List<String> nodeList = new ArrayList<String>(); //BFS visit list
//
//        Iterator leaves = net.getLeaves().iterator();
//        while(leaves.hasNext()){
//            NetNode temp = (NetNode)leaves.next();
//            in2leaf.put(temp.getName(),temp.getName());
//            heights.put(temp.getName(),0.0);
//        }
//        LinkedList q = new LinkedList<NetNode>();
//        NetNode root = net.getRoot();
//        nodeList.add(root.getName());
//        dist2Root.put(root.getName(),0.0);
//        Iterator child = root.getChildren().iterator();
//
//        while(child.hasNext()) {
//            NetNode c =(NetNode) child.next();
//            q.offer(c);
//            dist2Root.put(c.getName(),c.getParentDistance(root));
//        }
//        int cur ,last;
//        NetNode current = new BniNetNode();
//        double treeHeight = 0.0;
//        while(!q.isEmpty()) {
//            cur = 0;
//            last = q.size();
//            current = (NetNode) q.poll();
//            if (!current.isLeaf()) {
//                Iterator cTemp = current.getChildren().iterator();
//                while (cTemp.hasNext()) {
//                    q.offer(cTemp.next());
//                }
//            }
//            nodeList.add(current.getName());
//            NetNode parent = (NetNode) current.getParents().iterator().next();
//            double height2Root = (Double) dist2Root.get(parent.getName()) + current.getParentDistance(parent);
//            //cur++;
//            if (dist2Root.get(current.getName()) != null) {
//                if (height2Root < (Double) dist2Root.get(current.getName()))
//                    dist2Root.put(current.getName(), height2Root);
//            } else
//                dist2Root.put(current.getName(), height2Root);
//
//            if (height2Root > treeHeight)
//                treeHeight = height2Root;
//        }
//        //Iterator leafIt = nodeList.iterator();
//        List visitList = new ArrayList<String>();
//        for(int i = nodeList.size(); i>1;i--) {
//            NetNode tempNode = net.findNode((String)nodeList.get(i-1));
//            NetNode tempP = (NetNode) tempNode.getParents().iterator().next();
//            Iterator it = tempP.getChildren().iterator();
//            NetNode sib = (NetNode) it.next();
//            if (sib.equals(tempNode))
//                sib = (NetNode) it.next();
//            if(visitList.contains(tempNode.getName())||visitList.contains(sib.getName()))
//                continue;
//            double h1 = (double)heights.get(tempNode.getName());
//            double h2 = (double)heights.get(sib.getName());
//            heights.put(tempP.getName(), (treeHeight-(Double)dist2Root.get(tempP.getName()))/2);
//            if(h1<h2){
//                in2leaf.put(tempP.getName(),in2leaf.get(tempNode.getName()));
//
//                msCMD += " -ej " + heights.get(tempP.getName()).toString() + " " + in2leaf.get(sib.getName()) + " " + in2leaf.get(tempNode.getName());
//
//            }
//            else {
//                in2leaf.put(tempP.getName(), in2leaf.get(sib.getName()));
//                //msCMD += " -ej " + heights.get(tempP.getName()).toString() + " " + in2leaf.get(tempNode.getName()) + " " + in2leaf.get(sib.getName());
//                msCMD += " -ej " + heights.get(tempP.getName()).toString() + " " + in2leaf.get(tempNode.getName()) + " " + in2leaf.get(sib.getName());
//            }
//            visitList.add(tempNode.getName());
//            visitList.add(sib.getName());
//            //msCMD += " -ej " + heights.get(tempP.getName()).toString() + " " + in2leaf.get(sib.getName()) + " " + in2leaf.get(tempNode.getName());
//
//        }
////        for(int i = nodeList.size(); i>1;i--){
////
////            if(visitList.contains(nodeList.get(i-1)))
////                continue;
////            else{
////                NetNode tempNode = net.findNode(nodeList.get(i-1));
////                NetNode tempP = (NetNode)tempNode.getParents().iterator().next();
////                Iterator it = tempP.getChildren().iterator();
////                NetNode sib = (NetNode) it.next();
////                if (sib.equals(tempNode))
////                    sib = (NetNode) it.next();
////                msCMD += " -ej " + heights.get(tempP.getName()) + " " + in2leaf.get(sib.getName()) + " " + in2leaf.get(tempNode.getName());
////            }
////
////        }
//        System.out.println("MSCMD : " + msCMD.toString());
//        return msCMD;
//    }
//

    //nodesName can be used if leafName != int number
    public static String generateMSCMD(Network net, int numGT, Map<String, Integer> nodesName, HashMap<String, Double> cHeights) {
        //List<Tree> gts =  new ArrayList<Tree>();
        //net without outGroup
        Networks.autoLabelNodes(net);
        String msCMD = " " + net.getLeafCount() + " " + numGT + " -T -I" + " " + net.getLeafCount();
        for (int i = 0; i < net.getLeafCount(); i++) {
            msCMD += " 1";
        }
        //System.out.println("MSCMD : " + msCMD.toString());

        Map dist2Root = new HashMap<String, Double>();
        Map<String, Double> heights = new HashMap<String, Double>();
        Map in2leaf = new HashMap<String,String>();  //name mapping
        List<String> nodeList = new ArrayList<String>(); //BFS visit list

        Iterator leaves = net.getLeaves().iterator();
        while(leaves.hasNext()){
            NetNode temp = (NetNode)leaves.next();
            in2leaf.put(temp.getName(),temp.getName());
            heights.put(temp.getName(),0.0);
        }
        in2leaf.put("O", String.valueOf(net.getLeafCount()));
        LinkedList q = new LinkedList<NetNode>();
        NetNode root = net.getRoot();
        nodeList.add(root.getName());
        dist2Root.put(root.getName(),0.0);
        Iterator child = root.getChildren().iterator();

        while(child.hasNext()) {
            NetNode c =(NetNode) child.next();
            q.offer(c);
            dist2Root.put(c.getName(),c.getParentDistance(root));
        }
        int cur ,last;
        NetNode current = new BniNetNode();
        double treeHeight = 0.0;
        while(!q.isEmpty()) {
            cur = 0;
            last = q.size();
            current = (NetNode) q.poll();
            if (!current.isLeaf()) {
                Iterator cTemp = current.getChildren().iterator();
                while (cTemp.hasNext()) {
                    q.offer(cTemp.next());
                }
            }
            nodeList.add(current.getName());
            NetNode parent = (NetNode) current.getParents().iterator().next();
            double height2Root = (Double) dist2Root.get(parent.getName()) + current.getParentDistance(parent);
            //cur++;

                dist2Root.put(current.getName(), height2Root);

            if (height2Root > treeHeight)
                treeHeight = height2Root;
        }
        //Iterator leafIt = nodeList.iterator();
        List visitList = new ArrayList<String>();
        for(int i = nodeList.size(); i>1;i--) {
            NetNode tempNode = net.findNode((String)nodeList.get(i-1));
            NetNode tempP = (NetNode) tempNode.getParents().iterator().next();
            Iterator it = tempP.getChildren().iterator();
            NetNode sib = (NetNode) it.next();
            if (sib.equals(tempNode))
                sib = (NetNode) it.next();
            if(visitList.contains(tempNode.getName())||visitList.contains(sib.getName()))
                continue;
            double h1 = (double)heights.get(tempNode.getName());
            double h2 = (double)heights.get(sib.getName());
            heights.put(tempP.getName(), (treeHeight-(Double)dist2Root.get(tempP.getName()))/2);
            if(h1<h2){
                in2leaf.put(tempP.getName(),in2leaf.get(tempNode.getName()));

                msCMD += " -ej " + heights.get(tempP.getName()).toString() + " " + in2leaf.get(sib.getName()) + " " + in2leaf.get(tempNode.getName());

            }
            else {
                in2leaf.put(tempP.getName(), in2leaf.get(sib.getName()));
                msCMD += " -ej " + heights.get(tempP.getName()).toString() + " " + in2leaf.get(tempNode.getName()) + " " + in2leaf.get(sib.getName());
            }
            visitList.add(tempNode.getName());
            visitList.add(sib.getName());
            }
        //System.out.println("MSCMD : " + msCMD.toString());
        return msCMD;
    }



    public static Map<String, Double> getHeight(Network net, NetNode node, Map<String, Double> times, Map<String, String> interToLeaf) {
        double height = 0.0;
        List hList = new ArrayList<Double>();
        if (times.get(node.getName()) != null)
            height = times.get(node.getName());
        else {
            Iterator c = node.getChildren().iterator();
            String leafName = "";
            while (c.hasNext()) {
                NetNode child = (NetNode) c.next();
                getHeight(net, child, times, interToLeaf);
                if (times.get(child.getName()) + child.getParentDistance(node) > height) {
                    height = times.get(child.getName()) + child.getParentDistance(node);
                    leafName = interToLeaf.get(child.getName());
                    if (child.getName().equals(""))
                        System.out.print("Wrong");
                }
            }
            times.put(node.getName(), height);
            interToLeaf.put(node.getName(), leafName);
        }
        return times;
    }

    public static List<Tree> generateGTS(Network network, int numGTs, String MSPath, Map<String, String> MSName2AlleleName, HashMap<String, Double> heights, boolean ifOutGroup) {
        List<Tree> gts = new ArrayList<Tree>();
        Map<String, Integer> nodesName = new HashMap<String, Integer>();
        String cmd = generateMSCMD(network, numGTs, nodesName, heights);
        //Map<String, String> MSName2AlleleName = new HashMap<String, String>();
        String ogName = String.valueOf(network.getLeafCount());
        if(MSName2AlleleName == null){
            MSName2AlleleName = new HashMap<String, String>();
            Iterator leafIt = network.getLeaves().iterator();
            while(leafIt.hasNext()){
                String name = ((NetNode)leafIt.next()).getName();
                MSName2AlleleName.put(name,name);
            }
        }
        Iterator it = nodesName.keySet().iterator();
        while (it.hasNext()) {
            String name = (String) it.next();
            int num = nodesName.get(name);
            MSName2AlleleName.put(String.valueOf(num), name);
        }
        try {
            Process proc = Runtime.getRuntime().exec(MSPath + cmd, null, null);
            BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            String line;
            while ((line = stdout.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("//")) {
                    line = stdout.readLine();
                    String newLine = line;
                    if (ifOutGroup) {
                        int ogIndex = line.indexOf(ogName + ":");
                        newLine = line.substring(0, ogIndex) + "O:" + line.substring(ogIndex + ogName.length() + 1);

                    }
                    NewickReader nr = new NewickReader(new StringReader(newLine));
                    Tree newtr = nr.readTree();


                    for (TNode node : newtr.postTraverse()) {
                        if (node.isLeaf()) {
                            ((STINode) node).setName("MSTemp" + MSName2AlleleName.get(node.getName()));
                        }
                        //MS returns network with branch lengths in 4*N*mu, in Phylonet it is 2*N*mu
                        node.setParentDistance(node.getParentDistance() * 2);
                    }
                    for (TNode node : newtr.postTraverse()) {
                        if (node.isLeaf()) {
                            ((STINode) node).setName(node.getName().substring(6));
                        }
                    }
                    gts.add(newtr);
                }
            }

            proc.waitFor();
            stdout.close();
        } catch (Exception e) {
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return gts;
    }

    public static List<Tree> getMSTreeTopo(List<Tree> msGTS) {
        List<Tree> topos = new ArrayList<Tree>();
        for(int i=0;i< msGTS.size();i++){
            Tree temp = Trees.readTree(msGTS.get(i).toNewick());
            //
            msGTS.set(i,temp);
        }
        return  msGTS;
    }
}

//
//
// private static void processNexusFiles(int numGT) {
//        for(double scale : _scale



//    //Simulate by MS
//    protected List<Tree> simGTByMS(Network stTree, Integer numGTs){
//        SimGTInNetworkByMS simulator = new SimGTInNetworkByMS();
//        List<Tree> gts = new ArrayList<Tree>();
//        for (Tree tr : simulator.generateGTs(stTree, null, numGTs, _msdir)) {
//                gts.add(tr);
//            }
//        try {
//            PrintWriter out = new PrintWriter(numGTs + "GT.tree");
//            Iterator<Tree> it = gts.iterator();
//            while(it.hasNext()){
//                out.write(it.next().toString());
//            }
//            out.flush();
//            out.close();
//        } catch (IOException io) {
//            io.printStackTrace();
//        }
//        return gts;
//    }


//    //Simulate by Seq-gen
//    protected List<String> simSByGT(List gtTrees, Integer seqLength, int loci, double mutationRate) throws IOException{
//        List<String> seqs = new ArrayList<String>();
//        //seq-gen -mHKY -l 40 -s .2 <treefile >seqfile
//        String treeFile = gtTrees.size() + "GT.tree";
//        String seqFile = seqLength + "bp" + gtTrees.size() + "GT.txt";
//        Process proc = Runtime.getRuntime().exec(_seqdir+_seqGenCommand + "-l" + seqLength + "-s" + Double.toString(mutationRate) + "<" + treeFile + ">" + seqFile,null,null);
//        BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
//        String line;
//        while((line=stdout.readLine())!=null){
//            seqs.add(line);
//            }
//        return seqs;
//    }


//    public void generateFiles(String name) {
//        try {
//            PrintWriter out = new PrintWriter("./" + name + "/" + name + ".txt");
//            out.close();
//        } catch (IOException ioe) {
//            ioe.printStackTrace();
//        }
//    }

