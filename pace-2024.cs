using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pace_2024
{
   class Program
   {
      static void Main(string[] args)
      {
         PACE2024Solver solver = new PACE2024Solver();
         if (Environment.OSVersion.ToString().Contains("Windows"))
         {
            if (args.Length == 0)
            {
               int fnumber = 10;
               solver.RunSolver(fnumber.ToString() + ".gr");
            }
            else
            {
               solver.RunSolver(args[0], Convert.ToInt32(args[1]));
            }
         }
         else { solver.RunSolver(""); }
      }


      class PACE2024Solver
      {
         // control parameters:
         // the internal check terminates the program once a "runtime >= (MaxRuntime - 1)" is reached
         private const int MaxRuntime = 300; // TODO on optil.io: set this to around 294, 300 otherwise

         // data, from the actual instance and additionally generated:
         private int ASetSize, BSetSize, NoOfEdges;
         private Node[] AllNodes;
         private List<Edge> AllEdges;
         private List<Node> IsolatedBNodes;
         private List<Node> ActiveBNodes; // the 'active' B-nodes, i.e., the non-isolated ones
         private int NoOfActiveBNodes;
         private List<int>[] Neighbors; // for each node in the set B, the neighbors in set A connected by the edges
         private int[][] NeighborsARR; // the same as the list-datastructure, but in an array, which is faster to handle
         private bool MatricesAvailable;
         private int[,] XMatrix; // the number of crossings if node i precedes j,  \forall i,j \in B
         private int[,] ReversalCosts;
         private bool[,] Precedes;
         private long LB = 0;
         private List<IndependentInterval> IIs = new List<IndependentInterval>();

         // additional parameters:
         private const bool DisplayMessages = false; // TODO: this must be set to constant, = 'false'
         private const int DefaultRandomNumberSeed = 12; // default is 12
         private const int RevCostThreshold = 200; // good choice: 200 (experimentally verified)
         private const int ShiftCheck = 1300; // good choice: 1300 (experimentally verified)
         private const int MaxShift = 3900; // good choice: 3900, i.e., 3xShiftCheck (experimentally verified)
         private const double snmEQAcc = 0.5; // good choice = 0.5 (experimentally verified)
         private const double bmEQAcc = 0.1; // good choice = 0.1 (experimentally verified)
         private const double PermutationInversionPercentage = 0.2; // good choice = 0.2 (experimentally verified)
                                                                    //private const string ThePath = @"F:\PACE2024\"; // only used for the experiments on my machine
         private const int MatrixBound = 16384; // max size of the crossings matrix, i.e.: 16384 x 16384
                                                // other parameters (automatically set)
         private bool RunOnLinuxServer; // true if the code is run on Linux, false otherwise (= MS Windows)
         Random r; // random number, seed is given above
         DateTime OptimizationStart; // record the time (UTC) when optimization starts
         private long LoopCounter = 0; // counts the loops, required for initiating termination-checks
         private int PerturbationCounter = 0; // counts the perturbations (for statistics)
         private string TheFilename, TheSeed;


         public void RunSolver(string fn, int seed = DefaultRandomNumberSeed)
         {
            OptimizationStart = DateTime.UtcNow;
            System.Threading.Thread.CurrentThread.Priority = System.Threading.ThreadPriority.AboveNormal;
            TheFilename = fn;
            TheSeed = seed.ToString();

            if (Environment.OSVersion.ToString().Contains("Windows")) { RunOnLinuxServer = false; }
            else { RunOnLinuxServer = true; }

            r = new Random(seed);

            if (DisplayMessages)
            {
               Console.WriteLine("entering RunSolver_ILS for instance " + fn + " with seed=" + seed);
               Console.WriteLine("Run start=" + OptimizationStart.ToString());
            }

            if (DisplayMessages) { Console.WriteLine("Reading data..."); }
            ReadData(fn);
            if (DisplayMessages) { Console.WriteLine("Preprocessing..."); }
            Preprocessing();
            if (DisplayMessages) { Console.WriteLine("Computing independent intervals..."); }
            ComputeIndependentIntervals();

            if (MatricesAvailable)
            {
               if (DisplayMessages) { Console.WriteLine("Applying reductions..."); }
               Reductions();
            }
            if (DisplayMessages) { Console.WriteLine("Creating barycenter solution..."); }
            CreateBarycenterSolution();
            if (DisplayMessages) { Console.WriteLine("Computing bounds..."); }
            BoundsComputation();

            // initialize the elite solution
            foreach (IndependentInterval ii in IIs) { ii.UpdateElite(); }


            // check if the first solution is an optimum already
            if (GetEliteCosts() == LB)
            {
               if (DisplayMessages) { Console.WriteLine("barycenter solution is optimal"); }
               WriteEliteSolutionAndExitProgram(); return;
            }


            // a first run to the local optimum
            if (MatricesAvailable)
            {
               foreach (IndependentInterval ii in IIs)
               {
                  while (SingleNodeMove(ii)) { }
                  ii.UpdateElite();
               }
            }
            else
            {
               while (RunToLocalOptimum(IIs[0])) { };
               IIs[0].UpdateElite();
            }

            if (DisplayMessages) { Console.WriteLine("first local optimum=" + GetEliteCosts()); }
            // check if the first local optimum is a global optimum
            if (GetEliteCosts() == LB)
            {
               if (DisplayMessages) { Console.WriteLine("LB reached with first local optimum, LB = " + LB.ToString()); }
               WriteEliteSolutionAndExitProgram(); return;
            }


            // ILS until time runs out
            while (true)
            {
               if (MatricesAvailable)
               {
                  foreach (IndependentInterval ii in IIs)
                  {
                     if (ii.Gap() == 0) { continue; }
                     ii.ReturnToElite();
                     Perturbation(ii);
                     bool someImprovement = true;
                     while (someImprovement)
                     {
                        someImprovement = false;
                        while (SingleNodeMove(ii)) { someImprovement = true; }
                        while (BlockMove(ii, 2)) { someImprovement = true; }
                        while (BlockMove(ii, 3)) { someImprovement = true; }
                        while (BlockMove(ii, 4)) { someImprovement = true; }
                        while (BlockMove(ii, 5)) { someImprovement = true; }
                     }
                     bool u = ii.UpdateElite();
                     if (DisplayMessages && u)
                     {
                        Console.BackgroundColor = ConsoleColor.Blue;
                        Console.WriteLine(Environment.NewLine + "elite costs=" + GetEliteCosts().ToString() + "                 (update in ii" + ii.ID.ToString() + ")");
                        Console.BackgroundColor = ConsoleColor.Black;
                        Console.Title = GetEliteCosts().ToString("#,##0") + "   " + RecomputeEliteCosts().ToString("#,##0");
                     }
                  }
                  // check if the current elite is a global optimum
                  if (GetEliteCosts() == LB)
                  {
                     if (DisplayMessages) { Console.WriteLine("LB reached, LB = " + LB.ToString()); }
                     WriteEliteSolutionAndExitProgram(); return;
                  }
               }
               else
               {
                  IIs[0].ReturnToElite();
                  Perturbation(IIs[0]);
                  RunToLocalOptimum(IIs[0]);
                  IIs[0].UpdateElite();
               }
            }
         }


         private void ReadData(string fn)
         {
            AllEdges = new List<Edge>();
            System.IO.StreamReader sr = null;
            if (fn != "") { sr = new System.IO.StreamReader(@".\data\" + fn); }

            while (true)
            {
               string line = "";
               if (fn != "") { line = sr.ReadLine(); }
               else { line = Console.ReadLine(); }

               string[] splitted = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
               if (splitted[0] == "p")
               {
                  ASetSize = Convert.ToInt32(splitted[2]);
                  BSetSize = Convert.ToInt32(splitted[3]);
                  NoOfEdges = Convert.ToInt32(splitted[4]);
                  for (int i = 0; i < NoOfEdges; i++)
                  {
                     if (fn != "") { line = sr.ReadLine(); }
                     else { line = Console.ReadLine(); }
                     splitted = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                     AllEdges.Add(new Edge() { ANodeID = Convert.ToInt32(splitted[0]), BNodeID = Convert.ToInt32(splitted[1]) });
                  }
                  break;
               }
            }
            if (fn != "") { sr.Dispose(); }
         }
         private void Preprocessing()
         {
            // create node-objects
            AllNodes = new Node[ASetSize + BSetSize + 1];
            for (int i = 1; i < ASetSize + BSetSize + 1; i++)
            {
               Node n = new Node() { ProblemID = i, InternalID = -1 };
               AllNodes[i] = n;
            }
            // assign the edges to the nodes
            foreach (Edge e in AllEdges)
            {
               AllNodes[e.ANodeID].Neighbors.Add(e.BNodeID);
               AllNodes[e.BNodeID].Neighbors.Add(e.ANodeID);
            }

            // update the B-nodes
            IsolatedBNodes = new List<Node>();
            ActiveBNodes = new List<Node>();
            for (int i = 0; i < BSetSize; i++)
            {
               Node bnode = AllNodes[ASetSize + 1 + i];
               if (bnode.Neighbors.Count == 0)
               {
                  bnode.AV = int.MaxValue;
                  IsolatedBNodes.Add(bnode);
               }
               else
               {
                  bnode.AV = bnode.Neighbors.Average();
                  ActiveBNodes.Add(bnode);
               }
            }
            NoOfActiveBNodes = ActiveBNodes.Count();
            Neighbors = new List<int>[NoOfActiveBNodes];
            ActiveBNodes.Sort(CompareNodes); // sort the active B-nodes, then do the internal renumbering
            for (int i = 0; i < NoOfActiveBNodes; i++)
            {
               Node bnode = ActiveBNodes[i];
               bnode.InternalID = i;
               bnode.Neighbors.Sort();
               Neighbors[i] = bnode.Neighbors;
            }
            NeighborsARR = new int[BSetSize][];
            for (int i = 0; i < NoOfActiveBNodes; i++)
            {
               Node bnode = ActiveBNodes[i];
               int[] narr = new int[bnode.Neighbors.Count];
               NeighborsARR[bnode.InternalID] = narr;
               for (int pos = 0; pos < bnode.Neighbors.Count; pos++)
               { narr[pos] = bnode.Neighbors[pos]; }
            }

            // compute the crossing matrix
            MatricesAvailable = false;
            if (NoOfActiveBNodes > MatrixBound)
            {
               if (DisplayMessages) { Console.WriteLine("Matrix size above threshold."); }
               return;
            }
            try
            {
               XMatrix = new int[NoOfActiveBNodes, NoOfActiveBNodes];
               ReversalCosts = new int[NoOfActiveBNodes, NoOfActiveBNodes];
               Precedes = new bool[NoOfActiveBNodes, NoOfActiveBNodes];
            }
            catch
            {
               if (DisplayMessages) { Console.WriteLine("Unable to allocate matrices."); }
               return;
            }

            if (DisplayMessages) { Console.WriteLine("Matrices allocated...!"); }

            MatricesAvailable = true;
            int LeftNodeID, RightNodeID;
            int CrossCount, NeutralCount, ReverseCrossCount;
            LB = 0;

            for (LeftNodeID = 0; LeftNodeID < NoOfActiveBNodes - 1; LeftNodeID++)
            {
               int[] LeftNeighborList = NeighborsARR[LeftNodeID];
               for (RightNodeID = LeftNodeID + 1; RightNodeID < NoOfActiveBNodes; RightNodeID++)
               {
                  CrossCount = 0; // reset the count-variables
                  NeutralCount = 0;
                  int[] RightNeighborList = NeighborsARR[RightNodeID];
                  for (int pos1 = LeftNeighborList.Length - 1; pos1 > -1; pos1--)
                  {
                     int neighb1 = LeftNeighborList[pos1];
                     for (int pos2 = 0; pos2 < RightNeighborList.Length; pos2++)
                     {
                        int neighb2 = RightNeighborList[pos2];
                        if (neighb1 > neighb2) { CrossCount++; continue; }
                        else if (neighb1 == neighb2) { NeutralCount++; break; }
                        else { break; }
                     }
                  }

                  ReverseCrossCount = Neighbors[LeftNodeID].Count * Neighbors[RightNodeID].Count - CrossCount - NeutralCount;
                  XMatrix[LeftNodeID, RightNodeID] = CrossCount;
                  XMatrix[RightNodeID, LeftNodeID] = ReverseCrossCount;
                  ReversalCosts[LeftNodeID, RightNodeID] = ReverseCrossCount - CrossCount;
                  ReversalCosts[RightNodeID, LeftNodeID] = CrossCount - ReverseCrossCount;
                  LB += Math.Min(CrossCount, ReverseCrossCount);
               }
            }
         }
         private void ComputeIndependentIntervals()
         {
            //compute independent intervals
            if (!MatricesAvailable) // if there are no matrices, assign all B-nodes to a single interval
            {
               IndependentInterval oneInterval = new IndependentInterval() { Left = 1, Right = ASetSize, ID = 1 };
               IIs.Add(oneInterval);
               for (int i = 0; i < NoOfActiveBNodes; i++)
               {
                  Node bnode = ActiveBNodes[i]; //AllNodes[ASetSize + 1 + i];
                  bnode.IndependentIntervalID = 1;
                  oneInterval.Bnodes.Add(bnode.ProblemID);
               }
            }
            else
            {
               List<IndependentInterval> tmpcollection = new List<IndependentInterval>();
               int apos = 1;
               IndependentInterval CurrentInterval = new IndependentInterval() { Left = 1, Right = 1, ID = 1 };
               tmpcollection.Add(CurrentInterval);
               while (apos < ASetSize + 1)
               {
                  Node anode = AllNodes[apos];
                  //check if the rightmost position of the current internval has been reached
                  if (apos == CurrentInterval.Right)
                  {
                     //create the next interval
                     IndependentInterval NextInterval = new IndependentInterval() { Left = CurrentInterval.Right, Right = CurrentInterval.Right, ID = CurrentInterval.ID + 1 };
                     tmpcollection.Add(NextInterval);
                     foreach (int bnodeID in anode.Neighbors)
                     {
                        Node bnode = AllNodes[bnodeID];
                        if (bnode.Neighbors.Last() <= CurrentInterval.Right)
                        {
                           if (bnode.IndependentIntervalID < 0)
                           {
                              CurrentInterval.Bnodes.Add(bnodeID);
                              bnode.IndependentIntervalID = CurrentInterval.ID;
                           }
                        }
                        else
                        {
                           if (bnode.IndependentIntervalID < 0)
                           {
                              NextInterval.Bnodes.Add(bnodeID);
                              bnode.IndependentIntervalID = NextInterval.ID;
                           }
                           // check if the upcoming interval must be extended
                           if (bnode.Neighbors.Last() > NextInterval.Right)
                           {
                              NextInterval.Right = bnode.Neighbors.Last();
                           }
                        }
                     }
                     // continue with the next interval
                     CurrentInterval = NextInterval;
                  }
                  else // still in the interior of the interval
                  {
                     foreach (int bnodeID in anode.Neighbors)
                     {
                        Node bnode = AllNodes[bnodeID];
                        if (bnode.IndependentIntervalID < 0)
                        {
                           CurrentInterval.Bnodes.Add(bnodeID);
                           bnode.IndependentIntervalID = CurrentInterval.ID;
                        }
                        if (bnode.Neighbors.Last() > CurrentInterval.Right)
                        {
                           CurrentInterval.Right = bnode.Neighbors.Last();
                        }
                     }
                  }
                  apos++;
               }
               // consolidate the intervals
               foreach (IndependentInterval ai in tmpcollection)
               {
                  if (ai.Bnodes.Count > 0) { IIs.Add(ai); }
               }
            }
         }
         private int Reductions()
         {
            int ReductionCount = 0;
            // reduction rule RR1
            for (int LeftNodeID = 0; LeftNodeID < NoOfActiveBNodes - 1; LeftNodeID++)
            {
               for (int RightNodeID = LeftNodeID + 1; RightNodeID < NoOfActiveBNodes; RightNodeID++)
               {
                  if (XMatrix[LeftNodeID, RightNodeID] == 0 && XMatrix[RightNodeID, LeftNodeID] != 0) { Precedes[LeftNodeID, RightNodeID] = true; ReductionCount++; }
                  if (XMatrix[LeftNodeID, RightNodeID] != 0 && XMatrix[RightNodeID, LeftNodeID] == 0) { Precedes[RightNodeID, LeftNodeID] = true; ReductionCount++; }
               }
            }

            // reduction rule RR2: check for nodes with identical neighbors
            for (int pos1 = 0; pos1 < NoOfActiveBNodes - 1; pos1++)
            {
               Node node1 = ActiveBNodes[pos1];
               for (int pos2 = pos1 + 1; pos2 < NoOfActiveBNodes; pos2++)
               {
                  Node node2 = ActiveBNodes[pos2];
                  if (node1.AV != node2.AV) { break; }
                  if (!node1.NeighborsIdentical(node2)) { continue; }
                  //neighbors are identical for the two nodes, break the tie, make node1 precede node2
                  Precedes[node1.InternalID, node2.InternalID] = true;
                  ReductionCount++;
               }
            }
            return ReductionCount;
         }
         private void CreateBarycenterSolution()
         {
            // assign the barycenter-solution to the independent intervals
            // note that the active B-nodes are sorted before in the procedure 'Preprocessing'
            int currentpos = 0;
            foreach (IndependentInterval ie in IIs)
            {
               ie.permutation = new int[ie.Bnodes.Count];
               ie.RandomizedPositions = new int[ie.Bnodes.Count];
               ie.elite = new int[ie.Bnodes.Count];
               ie.eliteAdjacentCosts = long.MaxValue;
               for (int i = 0; i < ie.Bnodes.Count; i++)
               {
                  ie.permutation[i] = ActiveBNodes[currentpos].InternalID;
                  ie.RandomizedPositions[i] = i;
                  currentpos++;
               }
            }
         }
         private void BoundsComputation()
         {
            if (MatricesAvailable)
            {
               LB = 0;
               foreach (IndependentInterval ie in IIs)
               {
                  ie.UBstart = 0;
                  ie.LB = 0;
                  for (int pos1 = 0; pos1 < ie.Bnodes.Count - 1; pos1++)
                  {
                     int node1 = ie.permutation[pos1];
                     for (int pos2 = pos1 + 1; pos2 < ie.Bnodes.Count; pos2++)
                     {
                        int node2 = ie.permutation[pos2];
                        ie.UBstart += XMatrix[node1, node2];
                        ie.LB += Math.Min(XMatrix[node1, node2], XMatrix[node2, node1]);
                     }
                  }
                  ie.AdjacentCostsStart = 0;
                  for (int pos = 0; pos < ie.Bnodes.Count - 1; pos++)
                  {
                     int node1 = ie.permutation[pos];
                     int node2 = ie.permutation[pos + 1];
                     ie.AdjacentCostsStart += XMatrix[node1, node2];
                  }
                  ie.AdjacentCostsCurrent = ie.AdjacentCostsStart;
                  LB += ie.LB;
               }
            }
            else
            {
               foreach (IndependentInterval ie in IIs)
               {
                  ie.UBstart = Convert.ToInt64(NoOfEdges) * Convert.ToInt64(NoOfEdges);
                  ie.LB = 0;
                  ie.AdjacentCostsStart = 0;
                  for (int pos = 0; pos < ie.Bnodes.Count - 1; pos++)
                  {
                     int node1 = ie.permutation[pos];
                     int node2 = ie.permutation[pos + 1];
                     ie.AdjacentCostsStart += GetCrossings(node1, node2);
                  }
                  ie.AdjacentCostsCurrent = ie.AdjacentCostsStart;
               }
            }
         }


         private bool TimeoutCheck()
         {
            if (DisplayMessages) { Console.Write(" /c/ "); }
            if ((DateTime.UtcNow - OptimizationStart).TotalSeconds >= MaxRuntime - 1) { return true; }
            return false;
         }


         private void RandomizeRandomPositions(IndependentInterval ii)
         {
            for (int i = 0; i < ii.RandomizedPositions.Count(); i++)
            {
               int pos = r.Next(ii.RandomizedPositions.Count() - i);
               int numberAtpos = ii.RandomizedPositions[pos];
               ii.RandomizedPositions[pos] = ii.RandomizedPositions[ii.RandomizedPositions.Count() - 1 - i];
               ii.RandomizedPositions[ii.RandomizedPositions.Count() - 1 - i] = numberAtpos;
            }
         }


         private void Perturbation(IndependentInterval ii)
         {
            PerturbationCounter++;
            if (DisplayMessages) { Console.WriteLine(Environment.NewLine + "Perturbation"); }
            if (MatricesAvailable)
            {
               int MinInversionLength = 10;
               int MaxInversionLength = Convert.ToInt32(Math.Floor(Convert.ToDouble(ii.permutation.Length) * PermutationInversionPercentage));
               if (MinInversionLength > MaxInversionLength) { MinInversionLength = MaxInversionLength; }
               int InversionLength = MinInversionLength + r.Next(MaxInversionLength - MinInversionLength);
               int startpos = r.Next(ii.permutation.Length - InversionLength);
               long EvaluationChange = 0;
               for (int pos1 = startpos; pos1 < startpos + InversionLength - 1; pos1++)
               {
                  int node1 = ii.permutation[pos1];
                  for (int pos2 = pos1 + 1; pos2 < startpos + InversionLength; pos2++)
                  {
                     int node2 = ii.permutation[pos2];
                     EvaluationChange += ReversalCosts[node1, node2];
                  }
               }
               for (int i = 0; i < InversionLength / 2; i++)
               {
                  int pos1 = startpos + i;
                  int pos2 = startpos + InversionLength - 1 - i;
                  int node1 = ii.permutation[pos1];
                  int node2 = ii.permutation[pos2];
                  ii.permutation[pos1] = node2;
                  ii.permutation[pos2] = node1;
               }
               ii.AdjacentCostsCurrent += EvaluationChange;
            }
            else
            {
               int ReversalCount = ii.permutation.Length;
               if (DisplayMessages) { Console.WriteLine("perturbing in ii" + ii.ID.ToString() + " strength=" + ReversalCount.ToString()); }
               for (int i = 0; i < ReversalCount; i++)
               {
                  int pos1 = r.Next(ii.permutation.Length - 1);
                  int pos2 = pos1 + 1;
                  int node1 = ii.permutation[pos1];
                  int node2 = ii.permutation[pos2];
                  ii.AdjacentCostsCurrent += GetReversalDelta(node1, node2);
                  ii.permutation[pos1] = node2;
                  ii.permutation[pos2] = node1;
               }
            }
         }


         private bool BlockMove(IndependentInterval ii, int blocklength)
         {
            bool SomeImprovement = false;
            RandomizeRandomPositions(ii);
            for (int p = 0; p < ii.RandomizedPositions.Length; p++)
            {
               int blockstart = ii.RandomizedPositions[p];
               if (blockstart + blocklength > ii.permutation.Length) { continue; }
               LoopCounter++;
               if (LoopCounter % 10000 == 0 && TimeoutCheck()) { ii.UpdateElite(); WriteEliteSolutionAndExitProgram(); }
               int CumRevCosts = 0;
               int BestBlockShiftPos = -1;
               int BestBlockShiftCosts = int.MaxValue;
               for (int ShiftPos = blockstart + blocklength; ShiftPos < ii.permutation.Length; ShiftPos++)
               {
                  int NodeToTheRight = ii.permutation[ShiftPos];
                  bool StopMoving = false;
                  for (int i = 0; i < blocklength; i++)
                  {
                     int NodeInBlock = ii.permutation[blockstart + i];
                     if (Precedes[NodeInBlock, NodeToTheRight]) { StopMoving = true; break; }
                     CumRevCosts += ReversalCosts[NodeInBlock, NodeToTheRight];
                  }
                  if (StopMoving) { break; }
                  if (CumRevCosts > RevCostThreshold) { break; }
                  if (CumRevCosts < BestBlockShiftCosts) { BestBlockShiftCosts = CumRevCosts; BestBlockShiftPos = ShiftPos; }
               }
               CumRevCosts = 0;
               for (int ShiftPos = blockstart - 1; ShiftPos > -1; ShiftPos--)
               {
                  int NodeToTheLeft = ii.permutation[ShiftPos];
                  bool StopMoving = false;
                  for (int i = 0; i < blocklength; i++)
                  {
                     int NodeInBlock = ii.permutation[blockstart + i];
                     if (Precedes[NodeToTheLeft, NodeInBlock]) { StopMoving = true; break; }
                     CumRevCosts += ReversalCosts[NodeToTheLeft, NodeInBlock];
                  }
                  if (StopMoving) { break; }
                  if (CumRevCosts > RevCostThreshold) { break; }
                  if (CumRevCosts < BestBlockShiftCosts) { BestBlockShiftCosts = CumRevCosts; BestBlockShiftPos = ShiftPos; }
               }
               bool accept = false;
               if (BestBlockShiftCosts < 0) { accept = true; }
               else if (BestBlockShiftCosts == 0 && r.NextDouble() < bmEQAcc) { accept = true; }
               if (accept)
               {
                  if (BestBlockShiftPos > blockstart)
                  {
                     //move to the right
                     int[] block = new int[blocklength];
                     for (int i = 0; i < blocklength; i++) { block[i] = ii.permutation[blockstart + i]; }
                     int NoOfNodesThatMoveLeft = BestBlockShiftPos - blockstart - blocklength + 1;
                     int pos = blockstart;
                     for (int i = 0; i < NoOfNodesThatMoveLeft; i++)
                     {
                        ii.permutation[pos] = ii.permutation[pos + blocklength];
                        pos++;
                     }
                     for (int i = 0; i < blocklength; i++)
                     {
                        ii.permutation[pos] = block[i];
                        pos++;
                     }
                  }
                  else // move to the left
                  {
                     int[] block = new int[blocklength];
                     for (int i = 0; i < blocklength; i++) { block[i] = ii.permutation[blockstart + i]; }
                     int NoOfNodesThatMoveRight = blockstart - BestBlockShiftPos;
                     int pos = blockstart + blocklength - 1;
                     for (int i = 0; i < NoOfNodesThatMoveRight; i++)
                     {
                        ii.permutation[pos] = ii.permutation[pos - blocklength];
                        pos--;
                     }
                     pos = BestBlockShiftPos;
                     for (int i = 0; i < blocklength; i++)
                     {
                        ii.permutation[pos] = block[i];
                        pos++;
                     }
                  }
                  ii.AdjacentCostsCurrent += BestBlockShiftCosts;
                  if (BestBlockShiftCosts < 0) { SomeImprovement = true; }
                  p--;
               }
            }
            return SomeImprovement;
         }


         private long RecomputeEliteCosts()
         {
            if (!MatricesAvailable) { return long.MaxValue; }
            long c = 0;
            foreach (IndependentInterval ii in IIs)
            {
               for (int pos1 = 0; pos1 < ii.elite.Length - 1; pos1++)
               {
                  for (int pos2 = pos1 + 1; pos2 < ii.elite.Length; pos2++)
                  {
                     int nid1 = ii.elite[pos1];
                     int nid2 = ii.elite[pos2];
                     c += XMatrix[nid1, nid2];
                  }
               }
            }
            return c;
         }


         private bool SingleNodeMove(IndependentInterval ii)
         {
            bool SomeImprovement = false;
            RandomizeRandomPositions(ii);

            for (int p = 0; p < ii.RandomizedPositions.Length; p++)
            {
               int CurrentPos = ii.RandomizedPositions[p];
               LoopCounter++;
               if (LoopCounter % 10000 == 0 && TimeoutCheck()) { ii.UpdateElite(); WriteEliteSolutionAndExitProgram(); }

               int NodeToShift = ii.permutation[CurrentPos];
               int BestShiftPos = -1;
               int BestReversalCosts = int.MaxValue;
               // first look right
               int CumRevCosts = 0;
               for (int ShiftPos = CurrentPos + 1; ShiftPos < ii.permutation.Length; ShiftPos++)
               {
                  int NodeToTheRight = ii.permutation[ShiftPos];
                  if (Precedes[NodeToShift, NodeToTheRight]) { break; }
                  CumRevCosts += ReversalCosts[NodeToShift, NodeToTheRight];
                  if (CumRevCosts > RevCostThreshold) { break; }
                  if (CumRevCosts < BestReversalCosts) { BestReversalCosts = CumRevCosts; BestShiftPos = ShiftPos; }
               }
               // then look left
               CumRevCosts = 0;
               for (int ShiftPos = CurrentPos - 1; ShiftPos > -1; ShiftPos--)
               {
                  int NodeToTheLeft = ii.permutation[ShiftPos];
                  if (Precedes[NodeToTheLeft, NodeToShift]) { break; }
                  CumRevCosts += ReversalCosts[NodeToTheLeft, NodeToShift];
                  if (CumRevCosts > RevCostThreshold) { break; }
                  if (CumRevCosts < BestReversalCosts) { BestReversalCosts = CumRevCosts; BestShiftPos = ShiftPos; }
               }
               bool accept = false;
               if (BestReversalCosts < 0) { accept = true; }
               else if (BestReversalCosts == 0 && r.NextDouble() < snmEQAcc) { accept = true; }

               if (accept) // accept the shift move 
               {
                  if (BestShiftPos > CurrentPos) // make the shift to the right
                  {
                     for (int pos = CurrentPos; pos < BestShiftPos; pos++) { ii.permutation[pos] = ii.permutation[pos + 1]; }
                     ii.permutation[BestShiftPos] = NodeToShift;
                  }
                  else // make the shift to the left
                  {
                     for (int pos = CurrentPos; pos > BestShiftPos; pos--) { ii.permutation[pos] = ii.permutation[pos - 1]; }
                     ii.permutation[BestShiftPos] = NodeToShift;
                  }
                  // update costs
                  ii.AdjacentCostsCurrent += BestReversalCosts;
                  if (BestReversalCosts < 0)
                  {
                     //ImprovingMoves.Add("SNM " +loopcounter);
                     SomeImprovement = true;
                  }
                  p--;
               }
            }
            return SomeImprovement;
         }


         private long GetEliteCosts()
         {
            if (!MatricesAvailable && !RunOnLinuxServer)
            {
               long c = 0;
               foreach (IndependentInterval ii in IIs)
               {
                  c += ii.eliteAdjacentCosts - ii.AdjacentCostsStart;
               }
               if (TheFilename == "9.gr") { return 4300294 + c; }
               if (TheFilename == "10.gr") { return 18203862 + c; }
               if (TheFilename == "44.gr") { return 260871345686 + c; }
            }

            long tc = 0;
            foreach (IndependentInterval ii in IIs)
            {
               tc += (ii.UBstart - ii.AdjacentCostsStart + ii.eliteAdjacentCosts);
            }
            return tc;
         }


         private long GetTotalCosts()
         {

            if (!MatricesAvailable && !RunOnLinuxServer)
            {
               long c = 0;
               foreach (IndependentInterval ii in IIs)
               {
                  c += ii.AdjacentCostsCurrent - ii.AdjacentCostsStart;
               }
               if (TheFilename == "9.gr") { return 4300294 + c; }
               if (TheFilename == "10.gr") { return 18203862 + c; }
               if (TheFilename == "44.gr") { return 260871345686 + c; }
            }

            long tc = 0;
            foreach (IndependentInterval ii in IIs)
            {
               tc += (ii.UBstart - ii.AdjacentCostsStart + ii.AdjacentCostsCurrent);
            }
            return tc;
         }


         private bool RunToLocalOptimum(IndependentInterval ii)
         {
            if (DisplayMessages) { Console.WriteLine("entering RunToLocalOptimum with MaxShift=" + MaxShift); }
            RandomizeRandomPositions(ii);
            bool someimprovement = false;
            foreach (int CurrentPos in ii.RandomizedPositions)
            {
               LoopCounter++;
               if (LoopCounter % 1000 == 0)
               {
                  if (TimeoutCheck()) { ii.UpdateElite(); WriteEliteSolutionAndExitProgram(); }
                  if (DisplayMessages) { Console.Title = ("current = " + GetTotalCosts().ToString() + "  " + DateTime.UtcNow.ToString()); }
               }

               int NodeToShift = ii.permutation[CurrentPos];
               int MaxRightPosL1 = Math.Min(CurrentPos + 1 + ShiftCheck, ii.permutation.Length);
               int MaxRightPosL2 = Math.Min(CurrentPos + 1 + MaxShift, ii.permutation.Length);
               int MinLeftPosL1 = Math.Max(CurrentPos - 1 - ShiftCheck, -1);
               int MinLeftPosL2 = Math.Max(CurrentPos - 1 - MaxShift, -1);
               int BestShiftPos = -1;
               int BestReversalCosts = 0;
               int CumRevCostsToRight = 0;
               // first look right
               for (int ShiftPos = CurrentPos + 1; ShiftPos < MaxRightPosL1; ShiftPos++)
               {
                  CumRevCostsToRight += GetReversalDelta(NodeToShift, ii.permutation[ShiftPos]);
                  if (CumRevCostsToRight < BestReversalCosts)
                  {
                     BestReversalCosts = CumRevCostsToRight; BestShiftPos = ShiftPos;
                  }
               }
               int CumRevCostsToLeft = 0;
               // then left
               for (int ShiftPos = CurrentPos - 1; ShiftPos > MinLeftPosL1; ShiftPos--)
               {
                  CumRevCostsToLeft += GetReversalDelta(ii.permutation[ShiftPos], NodeToShift);
                  if (CumRevCostsToLeft < BestReversalCosts)
                  {
                     BestReversalCosts = CumRevCostsToLeft; BestShiftPos = ShiftPos;
                  }
               }
               // check if there is an improvement
               if (BestReversalCosts >= 0) { continue; } // if not, continue looking for one
                                                         // if so, then continue looking in the direction in which the improvement is found
               bool costsUpdated = false;
               int checkcounter = 0;
               if (BestShiftPos > CurrentPos) // look right
               {
                  for (int ShiftPos = MaxRightPosL1; ShiftPos < MaxRightPosL2; ShiftPos++)
                  {
                     CumRevCostsToRight += GetReversalDelta(NodeToShift, ii.permutation[ShiftPos]);
                     if (CumRevCostsToRight < BestReversalCosts)
                     {
                        BestReversalCosts = CumRevCostsToRight; BestShiftPos = ShiftPos;
                        costsUpdated = true;
                     }
                     checkcounter++;
                     if (checkcounter % ShiftCheck == 0 && !costsUpdated) { break; }
                  }
               }
               else // look left
               {
                  for (int ShiftPos = MinLeftPosL1; ShiftPos > MinLeftPosL2; ShiftPos--)
                  {
                     CumRevCostsToLeft += GetReversalDelta(ii.permutation[ShiftPos], NodeToShift);
                     if (CumRevCostsToLeft < BestReversalCosts)
                     {
                        BestReversalCosts = CumRevCostsToLeft; BestShiftPos = ShiftPos;
                        costsUpdated = true;
                     }
                     if (checkcounter % ShiftCheck == 0 && !costsUpdated) { break; }
                  }
               }

               if (BestReversalCosts < 0) // accept the shift move 
               {
                  if (BestShiftPos > CurrentPos) // make the shift to the right
                  {
                     for (int pos = CurrentPos; pos < BestShiftPos; pos++) { ii.permutation[pos] = ii.permutation[pos + 1]; }
                     ii.permutation[BestShiftPos] = NodeToShift;
                  }
                  else // make the shift to the left
                  {
                     for (int pos = CurrentPos; pos > BestShiftPos; pos--) { ii.permutation[pos] = ii.permutation[pos - 1]; }
                     ii.permutation[BestShiftPos] = NodeToShift;
                  }
                  // update costs
                  ii.AdjacentCostsCurrent += BestReversalCosts;
                  someimprovement = true;
               }
            }
            return someimprovement;
         }


         private void WriteEliteSolutionAndExitProgram()
         {
            StringBuilder sb = new StringBuilder();
            foreach (Node n in IsolatedBNodes) { sb.AppendLine(n.ProblemID.ToString()); }
            for (int iipos = 0; iipos < IIs.Count; iipos++)
            {
               IndependentInterval ii = IIs[iipos];
               for (int npos = 0; npos < ii.elite.Length; npos++)
               {
                  int internalID = ii.elite[npos];
                  int BnodeID = ActiveBNodes[internalID].ProblemID;
                  if (iipos == IIs.Count - 1 && npos == ii.elite.Length - 1)
                  {
                     sb.Append(BnodeID.ToString());
                  }
                  else
                  {
                     sb.AppendLine(BnodeID.ToString());
                  }
               }
            }
            if (RunOnLinuxServer)
            {
               Console.WriteLine(sb.ToString());
            }
            else
            {
               double actualruntime = (DateTime.UtcNow - OptimizationStart).TotalSeconds;
               if (DisplayMessages)
               {
                  Console.WriteLine("terminating with elite = " + GetEliteCosts());
                  Console.WriteLine("Runtime = " + actualruntime.ToString());
               }
               // write the elite solution
               System.IO.StreamWriter sw = new System.IO.StreamWriter(@".\results\" + TheFilename + "." + TheSeed + "." + GetEliteCosts().ToString() + ".sol", false);
               sw.Write(sb.ToString());
               sw.Dispose();
               bool resultswritten = false;
               while (!resultswritten)
               {
                  try
                  {
                     System.IO.StreamWriter sr = new System.IO.StreamWriter(@".\results\results.txt", true);
                     sr.WriteLine(TheFilename + "," + TheSeed + "," + actualruntime.ToString() + "," + GetEliteCosts().ToString());
                     sr.Dispose();
                     resultswritten = true;
                  }
                  catch { System.Threading.Thread.Sleep(1000); }
               }
            }
            System.Environment.Exit(0);
         }


         private static int CompareNodes(Node node1, Node node2)
         {
            if (node1.AV > node2.AV) { return 1; }
            if (node1.AV < node2.AV) { return -1; }
            return 0;
         }



         private int GetReversalDelta(int LeftNodeID, int RightNodeID)
         {
            int CrossCount = 0;
            int NeutralCount = 0;
            int[] LeftNeighborList = NeighborsARR[LeftNodeID];
            int[] RightNeighborList = NeighborsARR[RightNodeID];
            for (int pos1 = LeftNeighborList.Length - 1; pos1 > -1; pos1--)
            {
               int neighb1 = LeftNeighborList[pos1];
               for (int pos2 = 0; pos2 < RightNeighborList.Length; pos2++)
               {
                  int neighb2 = RightNeighborList[pos2];
                  if (neighb1 < neighb2) { break; }
                  else if (neighb1 == neighb2) { NeutralCount++; break; }
                  else { CrossCount++; }
               }
            }
            return LeftNeighborList.Length * RightNeighborList.Length - CrossCount - CrossCount - NeutralCount; ;
         }

         private int GetCrossings(int LeftNodeID, int RightNodeID)
         {
            int count = 0;
            int[] LeftNeighborList = NeighborsARR[LeftNodeID];
            int[] RightNeighborList = NeighborsARR[RightNodeID];
            for (int pos1 = LeftNeighborList.Length - 1; pos1 > -1; pos1--)
            {
               int neighb1 = LeftNeighborList[pos1];
               for (int pos2 = 0; pos2 < RightNeighborList.Length; pos2++)
               {
                  int neighb2 = RightNeighborList[pos2];
                  if (neighb1 <= neighb2) { break; }
                  count++;
               }
            }
            return count;
         }
      }


      class IndependentInterval
      {
         public int ID;
         public List<int> Bnodes = new List<int>();
         public int Left, Right;
         public int[] permutation;
         public long UBstart;
         public long AdjacentCostsStart;
         public long LB;
         public long AdjacentCostsCurrent;
         public int[] elite;
         public long eliteAdjacentCosts = long.MaxValue;
         public int[] RandomizedPositions;


         public long Gap()
         {
            return UBstart - AdjacentCostsStart + eliteAdjacentCosts - LB;
         }
         public bool UpdateElite()
         {
            bool TrueUpdate = false;
            if (AdjacentCostsCurrent <= eliteAdjacentCosts)
            {
               if (AdjacentCostsCurrent < eliteAdjacentCosts) { TrueUpdate = true; }
               permutation.CopyTo(elite, 0);
               eliteAdjacentCosts = AdjacentCostsCurrent;
            }
            return TrueUpdate;
         }
         public void ReturnToElite()
         {
            AdjacentCostsCurrent = eliteAdjacentCosts;
            elite.CopyTo(permutation, 0);
         }
      }


      class Edge
      {
         public int ANodeID, BNodeID;
      }
      class Node
      {
         public int InternalID;
         public int ProblemID;
         public double AV;
         public List<int> Neighbors = new List<int>();
         public int IndependentIntervalID = -1;
         public string GetNeighborString()
         {
            string s = "( ";
            for (int pos = 0; pos < Neighbors.Count - 1; pos++)
            {
               s += Neighbors[pos].ToString() + ", ";
            }
            s += Neighbors[Neighbors.Count - 1].ToString() + " )";
            return s;
         }
         public bool NeighborsIdentical(Node OtherNode)
         {
            if (AV != OtherNode.AV) { return false; }
            if (Neighbors.Count != OtherNode.Neighbors.Count) { return false; }
            for (int pos = 0; pos < Neighbors.Count(); pos++)
            {
               if (Neighbors[pos] != OtherNode.Neighbors[pos]) { return false; }
            }
            return true;
         }
      }

   }


}
