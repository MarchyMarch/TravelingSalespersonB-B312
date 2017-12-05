using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;


namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            string[] results = new string[3];

            int remainingCities = Cities.Length;
            int numOfSolutions = 0;
            int statesCreated = 0;
            int statesNotExplored = 0;

            DateTime start = DateTime.Now;
            DateTime end = start.AddSeconds(time_limit / 1000);

            BBState state = createState();
            statesCreated++;
            state.setPriority(calcKey(remainingCities - 1, state.getLowerBound()));

            double bssfBound = createGreedyBSSF();




            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }


        #region BBStateClass
        //***************************************************************************************************
        //****************************** State Class for Branch and Bound ***********************************
        //***************************************************************************************************

        /// <summary>
        /// This class tracks the state for each node along the path in the branch and bound algorithm
        /// </summary>
        public class BBState
        {
            private ArrayList path;
            private double lowerBound;
            private double priority;
            private double[,] costMatrix;

            // Constructor
            public BBState(ref ArrayList path, ref double lowerBound, ref double[,] costMatrix)
            {
                this.path = path;
                this.lowerBound = lowerBound;
                this.priority = double.MaxValue;
                this.costMatrix = costMatrix;
            }

            // Modifiers
            public void addToPath(City city)
            {
                path.Add(city);
            }

            public void setPriority(double priority)
            {
                this.priority = priority;
            }

            public void setLowerBound(double lowerBound)
            {
                this.lowerBound = lowerBound;
            }

            // Getters
            public ArrayList getPath()
            {
                return path;
            }

            public double getLowerBound()
            {
                return lowerBound;
            }

            public double getPriority()
            {
                return priority;
            }

            public double[,] getCostMatrix()
            {
                return costMatrix;
            }
        }
        #endregion

        #region AdditionalMethods
        //***************************************************************************************************
        //******************************** Additional Methods ***********************************************
        //***************************************************************************************************

        double calcKey(int remainingCities, double lowerBound)
        {
            if(remainingCities < 1)
            {
                return lowerBound;
            }
            else
            {
                double result = (lowerBound / (Cities.Length - remainingCities));
                return result;
            }
        }

        double createGreedyBSSF()
        {
            Route = new ArrayList();
            Route.Add(Cities[0]);
            int curCityIndex = 0;

            while(Route.Count < Cities.Length)
            {
                double minValue = double.MaxValue;
                int minIndex = 0;

                for(int i = 0; i < Cities.Length; i++)
                {
                    if(curCityIndex != i)
                    {
                        if(!Route.Contains(Cities[i]))
                        {
                            double compareValue = Cities[curCityIndex].costToGetTo(Cities[i]);

                            if(compareValue < minValue)
                            {
                                if(Route.Count == Cities.Length-1 && Cities[i].costToGetTo(Cities[0]) == double.MaxValue)
                                {
                                    continue;
                                }

                                minValue = compareValue;
                                minIndex = i;
                            }
                        }
                    }
                }

                curCityIndex = minIndex;
                Route.Add(Cities[curCityIndex]);
            }

            bssf = new TSPSolution(Route);
            return bssf.costOfRoute();
        }

        void setUpMatrix(ref double[,] costMatrix, int parentIndex, int childIndex, ref double lowerBound)
        {
            if(costMatrix[parentIndex, childIndex] != double.MaxValue)
            {
                lowerBound += costMatrix[parentIndex, childIndex];
            }

            for(int row = 0; row < Cities.Length; row++)
            {
                costMatrix[row, childIndex] = double.MaxValue;
            }

            for(int col = 0; col < Cities.Length; col++)
            {
                costMatrix[parentIndex, col] = double.MaxValue;
            }

            costMatrix[parentIndex, childIndex] = double.MaxValue;
        }

        double reduceMatrix(ref double[,] costMatrix)
        {
            double lowerBound = 0;

            for(int row = 0; row < Cities.Length; row++)
            {
                double minValue = double.MaxValue;
                for(int col = 0; col < Cities.Length; col++)
                {
                    if(costMatrix[row,col] < minValue)
                    {
                        minValue = costMatrix[row, col];
                    }
                }

                if(minValue != 0 && minValue != double.MaxValue)
                {
                    lowerBound += minValue;

                    for(int col = 0; col < Cities.Length; col++)
                    {
                        if(costMatrix[row,col] != double.MaxValue)
                        {
                            costMatrix[row, col] -= minValue;
                        }
                    }
                }
            }

            return lowerBound;
        }

        BBState createState()
        {
            double[,] initialCostMatrix = new double[Cities.Length, Cities.Length];

            for(int i = 0; i < Cities.Length; i++)
            {
                for(int j = 0; j < Cities.Length; j++)
                {
                    if(i == j)
                    {
                        initialCostMatrix[i, j] = double.MaxValue;
                    }
                    else
                    {
                        initialCostMatrix[i, j] = Cities[i].costToGetTo(Cities[j]);
                    }
                }
            }

            ArrayList path = new ArrayList();
            path.Add(Cities[0]);

            double lowerBound = reduceMatrix(ref initialCostMatrix);

            BBState state = new BBState(ref path, ref lowerBound, ref initialCostMatrix);

            return state;
        }
        #endregion

        #region PriorityQueue
        //***************************************************************************************************
        //************************************ Priority Queue ***********************************************
        //***************************************************************************************************

        public sealed class PriorityQueue
        {
            private int capacity;
            private int count;
            private int maxNumber;
            private BBState[] states;

            public PriorityQueue() { }

            public bool isEmpty()
            {
                return count == 0;
            }

            public int getCount()
            {
                return count;
            }

            public double getMinLowerBound()
            {
                return states[0].getLowerBound();
            }

            public int getMaxNumber()
            {
                return maxNumber;
            }

            public void makeQueue(int numNodes)
            {
                states = new BBState[1000000];
                capacity = numNodes;
                count = 0;
                maxNumber = 0;
            }

            public BBState deleteMin()
            {
                BBState minState = states[1];

                states[1] = states[count];
                count--;

                int whileItterator = 1;

                while(whileItterator <= count)
                {
                    int leftChildIndex = whileItterator * 2;
                    if (leftChildIndex > count) break;

                    if(leftChildIndex+1 <= count && states[leftChildIndex +1].getPriority() < states[leftChildIndex].getPriority())
                    {
                        leftChildIndex++;
                    }

                    if(states[whileItterator].getPriority() > states[leftChildIndex].getPriority())
                    {
                        BBState tempState = states[leftChildIndex];
                        states[leftChildIndex] = states[whileItterator];
                        states[whileItterator] = tempState;
                    }

                    whileItterator = leftChildIndex;
                }

                return minState;
            }

            public void insert(BBState state)
            {
                count++;
                states[count] = state;
                if (count > maxNumber) maxNumber = count;

                int whileIterrator = count;
                while(whileIterrator > 1 && states[whileIterrator/2].getPriority() > states[whileIterrator].getPriority())
                {
                    BBState tempState = states[whileIterrator / 2];
                    states[whileIterrator / 2] = states[whileIterrator];
                    states[whileIterrator] = tempState;

                    whileIterrator = whileIterrator / 2;
                }
            }
        }

        #endregion

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];

            // TODO: Add your implementation for a greedy solver here.

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }

        public string[] fancySolveProblem()
        {
            string[] results = new string[3];

            // TODO: Add your implementation for your advanced solver here.

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }
        #endregion
    }

}
