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

        #region B&Bsolution
        //***************************************************************************************************
        //*********************************** Branch and Bound **********************************************
        //***************************************************************************************************

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// Time Complexity: Because of the while loop, which is the meat of the solution, the time complexity
        /// is O(2^N * N^2) 
        /// Space Complexity: Like the time complexity the main memory allocation happens in the while loop 
        /// and also has a complexity of O(2^N * N^2)
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, 
        /// time spent to find solution, number of solutions found during search 
        /// (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            string[] results = new string[3];

            // numbers we are looking for and help with solving the problem
            int remainingCities = Cities.Length;
            int numOfSolutions = 0;
            int statesCreated = 0;
            int statesNotExplored = 0;

            // the start of when the algorithm was called and the end being the added time limit
            DateTime start = DateTime.Now;
            DateTime end = start.AddSeconds(time_limit / 1000);

            // create the initial state at the start city then set its priority
            // the createState() method makes this part a time complexity of O(N^2)
            BBState state = createState();
            statesCreated++;
            state.setPriority(calcKey(remainingCities - 1, state.getLowerBound()));

            // creat the bssf with a greedy approach
            // Making a greedy BSSF gives this part a time complexity of O(N^2)
            // as well as a space complexity of O(N)
            double bssfBound = createGreedyBSSF();

            // makes a new priority queue and adds the start city
            // this block only allocates space in makeQueue which is O(1,000,000)
            PriorityQueue queue = new PriorityQueue();
            queue.makeQueue(Cities.Length);
            queue.insert(state);

            // This is the implementation of the branch and bound and finds the optimal solution
            // unless the time limit is met.  Here time complexity is terrible becasue it is 
            // O(2^N * N^2).  it itterates through the while loop O(2^N) and for each state it
            // reduces the matrix with a time complexity of O(N^2). For space complexity each state 
            // makes an N by N matrix giving it a space comlexity of O(N^2) for each matrix and it
            // also makes 2^N states so its overall space complexity is O(2^N * N^2)
            while(!queue.isEmpty() && DateTime.Now < end && queue.getMinLowerBound() != bssfBound)
            {
                // get the next city
                BBState curState = queue.deleteMin();

                // check to see if it has a better bssf
                if(curState.getLowerBound() < bssfBound)
                {
                    // branch from the state to explore its children
                    for(int i = 0; i < Cities.Length; i++)
                    {
                        if (DateTime.Now >= end) break;

                        // if we have visited the city we skip it
                        if (curState.getPath().Contains(Cities[i])) continue;

                        // start to make a new state
                        double[,] curCostMatrix = curState.getCostMatrix();
                        double[,] newCostMatrix = new double[Cities.Length, Cities.Length];

                        // fill the new matrix
                        // this is O(N^2)
                        for(int j = 0; j < Cities.Length; j++)
                        {
                            for(int k = 0; k < Cities.Length; k++)
                            {
                                newCostMatrix[j, k] = curCostMatrix[j, k];
                            }
                        }
                        City lastCityInCurState = (City)curState.getPath()[curState.getPath().Count - 1];
                        double curLowerBound = curState.getLowerBound();
                        setUpMatrix(ref newCostMatrix, Array.IndexOf(Cities, lastCityInCurState), i, ref curLowerBound);
                        double newLowerBound = curLowerBound + reduceMatrix(ref newCostMatrix);

                        ArrayList curPath = curState.getPath();
                        ArrayList newPath = new ArrayList();

                        foreach(City city in curPath)
                        {
                            newPath.Add(city);
                        }
                        newPath.Add(Cities[i]);

                        // makes new state
                        BBState childState = new BBState(ref newPath, ref newLowerBound, ref newCostMatrix);
                        statesCreated++;

                        // don't explore child states that larger than the bssf
                        if(childState.getLowerBound() < bssfBound)
                        {
                            City startCity = (City)childState.getPath()[0];
                            City finalCity = (City)childState.getPath()[childState.getPath().Count - 1];
                            double costToReturn = finalCity.costToGetTo(startCity);

                            // check to see if it goes to the start city
                            if(childState.getPath().Count == Cities.Length && costToReturn != double.MaxValue)
                            {
                                childState.setLowerBound(childState.getLowerBound() + costToReturn);
                                bssf = new TSPSolution(childState.getPath());
                                bssfBound = bssf.costOfRoute();
                                numOfSolutions++;
                                statesNotExplored++;
                            }
                            else
                            {
                                // set the priority and add it to the queue
                                remainingCities = Cities.Length - childState.getPath().Count;
                                childState.setPriority(calcKey(remainingCities, childState.getLowerBound()));
                                queue.insert(childState);
                            }
                        }
                        else
                        {
                            statesNotExplored++;
                        }
                    }
                }
                curState = null;
            }

            // add the queue length in case the algorithm times out
            statesNotExplored += queue.getCount();

            // find how long it took to find the optimal route and set the span to seconds
            end = DateTime.Now;
            TimeSpan timeSpan = end - start;
            double spanSeconds = timeSpan.TotalSeconds;

            // convert all the important variables that need to be returned to strings and return them
            results[COST] = System.Convert.ToString(bssf.costOfRoute());
            results[TIME] = System.Convert.ToString(spanSeconds);
            results[COUNT] = System.Convert.ToString(numOfSolutions);

            return results;
        }
        #endregion

        #region BBStateClass
        //***************************************************************************************************
        //****************************** State Class for Branch and Bound ***********************************
        //***************************************************************************************************

        /// <summary>
        /// This class tracks the state for each node along the path in the branch and bound algorithm
        /// Time Complexity:  All methods are O(1)
        /// Space Complexity: O(N^2) because each new state has space allocated for the costMatrix which has
        /// a space of an N by N matrix.
        /// The data structures i used were an array list to keep track of the path, for the cost matrix I 
        /// used a 2d array of doubles, and the other private members i used a double for large cases. Most
        /// types are primitive because Dr. Ventura said they are faster.
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

        #region PriorityQueue
        //***************************************************************************************************
        //************************************ Priority Queue ***********************************************
        //***************************************************************************************************

        /// <summary>
        /// Modified version of the priority queue from lab 03.  It uses an array to keep track of the states
        /// and uses the bubbling method to sort the array by highest priority.  The priority key is dependent 
        /// on both the lower bound and the cities we have visitied.  This shows which cities are closer and
        /// where we are to visit next
        /// </summary>
        public sealed class PriorityQueue
        {
            private int capacity;
            private int count;
            private int maxNumber;
            private BBState[] states;

            public PriorityQueue() { }
            
            // Checks if queue is empty
            public bool isEmpty()
            {
                return count == 0;
            }

            // Getters
            public int getCount()
            {
                return count;
            }

            public double getMinLowerBound()
            {
                return states[1].getLowerBound();
            }

            public int getMaxNumber()
            {
                return maxNumber;
            }

            /// <summary>
            /// makes the initial queue.  It initializes the states array to a million then sets the capacity
            /// to the number of nodes, which is the number of cities.
            /// Time Complexity: O(1), there are no loops or itterations, just assignments
            /// Space Complexity: 1,000,000 because thats what is initialized for the queue, it is a constant
            /// </summary>
            /// <param name="numNodes">Equal to the number of cities to explore</param>
            public void makeQueue(int numNodes)
            {
                states = new BBState[1000000];
                capacity = numNodes;
                count = 0;
                maxNumber = 0;
            }

            /// <summary>
            /// Here we delete the minimum state and then adjust the tree accordingly, similar to lab 03
            /// Time Complexity: the time complexity for this is O(logN) where N is the height of the tree
            /// Space Complexity: O(1) because there are no new space allocations
            /// </summary>
            /// <returns>The minimum BBState/optimal city</returns>
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

            /// <summary>
            /// This inserts the new state at the bottom of the tree and then it bubbles up to the correct place
            /// Time Complexity: the worse case scenario is if the added state is the root of the tree so the 
            /// time complexity is O(logN) where N is the height of the tree
            /// Space Complexity: because there is now new allocation of space the complexity is O(1)
            /// </summary>
            /// <param name="state">The state to be added to the queue</param>
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

        #region AdditionalMethods
        //***************************************************************************************************
        //******************************** Additional Methods ***********************************************
        //***************************************************************************************************


        /// <summary>
        /// Makes a key to a state for use in the priority queue
        /// Time Complexity: O(1) because there are only if statements and mathmatical computations
        /// Space Complexity: O(1) there is nothing that is stored through this method
        /// </summary>
        /// <param name="remainingCities">Cities Remaining to be visited</param>
        /// <param name="lowerBound">The current low cost</param>
        /// <returns>The key to the set priority</returns>
        double calcKey(int remainingCities, double lowerBound)
        {
            if (remainingCities < 1)
            {
                return lowerBound;
            }
            else
            {
                double result = (lowerBound / (Cities.Length - remainingCities));
                return result;
            }
        }


        /// <summary>
        /// Helps to make an initial BSSF
        /// Time Complexity: This method itterates through all the cities and for each city it itterates
        /// through all the cities again, so it has a time complexity of O(N^2)
        /// Space Complexity: This method makes an array list of all the cities so the space complexity
        /// is O(N)
        /// </summary>
        /// <returns>Initial BSSF value</returns>
        double createGreedyBSSF()
        {
            Route = new ArrayList();
            Route.Add(Cities[0]);
            int curCityIndex = 0;

            while (Route.Count < Cities.Length)
            {
                double minValue = double.MaxValue;
                int minIndex = 0;

                for (int i = 0; i < Cities.Length; i++)
                {
                    if (curCityIndex != i)
                    {
                        if (!Route.Contains(Cities[i]))
                        {
                            double compareValue = Cities[curCityIndex].costToGetTo(Cities[i]);

                            if (compareValue < minValue)
                            {
                                if (Route.Count == Cities.Length - 1 && Cities[i].costToGetTo(Cities[0]) == double.MaxValue)
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

        /// <summary>
        /// Modifies the lowerBound and then modifies the rows and columns of the cost matrix
        /// Time Complexity: because we only modify the row and column the time complexity is
        /// O(N) where N is the length of the row or column
        /// Space Complexity: O(1) everything that is modified is done so by reference so
        /// there is no new allocation of memory
        /// </summary>
        /// <param name="costMatrix">the given cost matrix</param>
        /// <param name="parentIndex">the index of parent city</param>
        /// <param name="childIndex">the index of the next city</param>
        /// <param name="lowerBound">current lower bound</param>
        void setUpMatrix(ref double[,] costMatrix, int parentIndex, int childIndex, ref double lowerBound)
        {
            if (costMatrix[parentIndex, childIndex] != double.MaxValue)
            {
                lowerBound += costMatrix[parentIndex, childIndex];
            }

            for (int row = 0; row < Cities.Length; row++)
            {
                costMatrix[row, childIndex] = double.MaxValue;
            }

            for (int col = 0; col < Cities.Length; col++)
            {
                costMatrix[parentIndex, col] = double.MaxValue;
            }

            costMatrix[parentIndex, childIndex] = double.MaxValue;
        }

        /// <summary>
        /// reduces the cost matrix and calculate lower bound
        /// Time Complexity: It itterates through the cost matrix, size N by N, twice it has a time 
        /// complexity of O(N^2)
        /// Space Complexity: There is no new allocation of space and the modified matrix is passed
        /// in by reference so the complexity is O(1)
        /// </summary>
        /// <param name="costMatrix">the matrix to reduce</param>
        /// <returns>the modified lower bound</returns>
        double reduceMatrix(ref double[,] costMatrix)
        {
            double lowerBound = 0;

            for (int row = 0; row < Cities.Length; row++)
            {
                double minValue = double.MaxValue;
                for (int col = 0; col < Cities.Length; col++)
                {
                    if (costMatrix[row, col] < minValue)
                    {
                        minValue = costMatrix[row, col];
                    }
                }

                if (minValue != 0 && minValue != double.MaxValue)
                {
                    lowerBound += minValue;

                    for (int col = 0; col < Cities.Length; col++)
                    {
                        if (costMatrix[row, col] != double.MaxValue)
                        {
                            costMatrix[row, col] -= minValue;
                        }
                    }
                }
            }

            return lowerBound;
        }

        /// <summary>
        /// creates a list for the first city in the list
        /// Time Complexity: O(N^2) because it itterates through the cities twice
        /// Space Complexity: Because spaceis allocated for the matrix space complexity is O(N^2)
        /// </summary>
        /// <returns></returns>
        BBState createState()
        {
            double[,] initialCostMatrix = new double[Cities.Length, Cities.Length];

            for (int i = 0; i < Cities.Length; i++)
            {
                for (int j = 0; j < Cities.Length; j++)
                {
                    if (i == j)
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
