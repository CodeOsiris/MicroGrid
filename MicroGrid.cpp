#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

#ifdef _WIN32
#include <windows.h>
#define VERSION 0
#else
#define VERSION 1
#endif

//Basic Const
const int POP_SIZE = 100;
const int MAX_GENERATION = 12000;
const int NUM_OF_GENE = 240;
const int MAX_STABLE = 800;
const int MODGA_R = 50;
const double P_CROSSOVER = 0.7;
const double P_MUTATION = 0.02;
const double P_SUBSTITUTE = 0.05;
const string INPUT_FILE_NAME = "ga_data";
const string OUTPUT_FILE_NAME = "ga_log";

//Particular Const
const int NUM_OF_CITY = 5;
const int GENE_PER_POINT = 10;
const int CHECK_POINT = 24;
const double PRICE_PER_UNIT = 1.0;
const double COST_PER_UNIT = 0.3;
const double COST_PER_DIST = 0.1;
const double LOST_PER_LACK = 1.0;
const double GAIN_PER_MORE = 0.05;
const double ROLL_BOUND = 0.01;
const int MAX_LOOP_TIME = 15;
const int CURRENT[GENE_PER_POINT] = {0,0,0,0,1,1,1,2,2,3};
const int BORROW[GENE_PER_POINT] = {1,2,3,4,2,3,4,3,4,4};

//Basic Variable
int generation;
int best_record;
int stable;
int best_mem;
int worst_mem;

//Particular Variable
int map[NUM_OF_CITY][NUM_OF_CITY] = {0};
int random_sort[GENE_PER_POINT];
double usage[NUM_OF_CITY][CHECK_POINT] = {0};
double generate[NUM_OF_CITY][CHECK_POINT] = {0};
double actual[NUM_OF_CITY][CHECK_POINT] = {0};
double storage[NUM_OF_CITY][CHECK_POINT] = {0};
double differ[NUM_OF_CITY][CHECK_POINT] = {0};
bool able[NUM_OF_CITY][CHECK_POINT] = {false};

//SA mutation
double adjust[MAX_GENERATION + 1];

class Gene
{
public:
	double borrow_amount;
};

class Genetype
{
public:
	Genetype():fitness(0),rfitness(0),cfitness(0){}
	//Basic Property
	Gene gene[NUM_OF_GENE];
	double upper[NUM_OF_GENE];
	double lower[NUM_OF_GENE];
	double fitness;
	double rfitness;
	double cfitness;
    double punish;
	bool chosen;
};
Genetype population[POP_SIZE + 1];
Genetype newpopulation[POP_SIZE + 1];
Genetype history_best;

int maxcity(int,int);
void initialize(int);
double random_value(double,double);
void evaluate();
void keep_best();
void elitist();
void select();
void crossover();
void crossover(int,int);
int compete(int,int);
void mutate();
void rollback();
void swirl(int,Genetype &);
void substitute();
void substitute(int,int);
void report(int);
void trace(int);

int maxcity(int a,int b)
{
	double amax = 0;
	int v = -1,point,current_city,borrow_city;
	for (int i = 0;i < GENE_PER_POINT;i++)
	{
		point = b;
		current_city = CURRENT[i];
		borrow_city = BORROW[i];
		if (storage[borrow_city][point] > amax)
		{
			amax = storage[borrow_city][point];
			v = BORROW[i];
		}
	}
	return v;
}

void init_file()
{
	ifstream fin("ga_data0.txt");
	int i,j,point,current_city,borrow_city;
	stable = 0;
	double sum = 0;
	//读取环境变量
	for (i = 0;i < NUM_OF_CITY;i++)
	{
		for (j = 0;j < CHECK_POINT;j++)
			fin >> usage[i][j];
		for (j = 0;j < CHECK_POINT;j++)
			fin >> generate[i][j];
	}
	for (i = 0;i < NUM_OF_CITY;i++)
		for (j = 0;j < NUM_OF_CITY;j++)
			fin >> map[i][j];
	//根据特定方案初始化
	ifstream finf("ga_total.txt");
	int city;
	double borrow;
	for (int i = 0;i < 24;i++)
	{
		for (int j = 0;j < 5;j++)
		{
			finf >> city >> borrow;
			for (int k = 0;k < 10;k++)
			{
				if (CURRENT[k] == j && BORROW[k] == city)
					population[0].gene[i * 10 + k].borrow_amount = borrow;
				else if (CURRENT[k] == city && BORROW[k] == j)
					population[0].gene[i * 10 + k].borrow_amount = -borrow;
			}
		}
	}
	evaluate();
	cout << population[0].fitness;
	cin >> city;
	finf.close();
}

void initialize(int cnt)
{
	/*
	for (int i = 0;i < POP_SIZE;i++)
		for (int j = 0;j < NUM_OF_GENE;j++)
		{
			population[i].lower[j] = -1.28;
			population[i].upper[j] = 1.28;
			population[i].gene[j].borrow_amount = random_value(population[i].lower[j],population[i].upper[j]);
		}
	population[POP_SIZE] = history_best = population[0];
	*/
	ifstream fin(INPUT_FILE_NAME + to_string(cnt) + ".txt");
	int i,j,point,current_city,borrow_city;
    bool flag[GENE_PER_POINT] = {0};
    random_sort[0] = rand() % GENE_PER_POINT;
    flag[random_sort[0]] = true;
    for (i = 1;i < GENE_PER_POINT;i++)
    {
        j = rand() % GENE_PER_POINT;
        while (flag[(random_sort[i - 1] + j) % GENE_PER_POINT])
            j = (j + rand() % GENE_PER_POINT) % GENE_PER_POINT;
        random_sort[i] = (random_sort[i - 1] + j) % GENE_PER_POINT;
        flag[random_sort[i]] = true;
    }

	stable = 0;
	for (i = 0;i < NUM_OF_CITY;i++)
	{
		for (j = 0;j < CHECK_POINT;j++)
			fin >> usage[i][j];
		for (j = 0;j < CHECK_POINT;j++)
			fin >> generate[i][j];
	}
    for (i = 0;i < NUM_OF_CITY;i++)
    {
        double sum[CHECK_POINT] = {0};
        for (j = 0;j < CHECK_POINT;j++)
        {
            if (j > 0)
              sum[j] = sum[j - 1];
            sum[j] += generate[i][j] - usage[i][j];
        }
        for (j = CHECK_POINT - 1;j >= 0;j--)
        {
            if (sum[j] < 0)
            {
                for (j = j;j >= 0;j--)
                  able[i][j] = false;
            }
            else able[i][j] = true;
        }
    }
	for (i = 0;i < NUM_OF_CITY;i++)
		for (j = 0;j < NUM_OF_CITY;j++)
			fin >> map[i][j];
	for (i = 0;i < POP_SIZE;i++)
	{	
		for (j = 0;j < NUM_OF_GENE;j++)
		{
			point = j / GENE_PER_POINT;
			current_city = CURRENT[j % GENE_PER_POINT];
			borrow_city = BORROW[j % GENE_PER_POINT];
			//Initialize storage,borrow_city,borrow_amount
			population[i].lower[j] = min(0.0,usage[current_city][point] - generate[current_city][point]);
			population[i].upper[j] = max(0.0,generate[borrow_city][point] - usage[borrow_city][point]);
            if (able[current_city][point] && able[borrow_city][point])
    			population[i].gene[j].borrow_amount = random_value(population[i].lower[j],population[i].upper[j]);
            else if (able[borrow_city][point])
              population[i].gene[j].borrow_amount = random_value(0,population[i].upper[j]);
            else if (able[current_city][point])
              population[i].gene[j].borrow_amount = random_value(population[i].lower[j],0);
            else population[i].gene[j].borrow_amount = 0;
		}
	}
	population[POP_SIZE] = population[0];
	fin.close();
}

double random_value(double low,double high)
{
	if (low == high)
		return low;
	if (low > high)
		swap(low,high);
	double value;
	value = (rand() % 1000) * (high - low) / 1000.0 + low;
	return value;
}

void evaluate()
{
	int i,j,k,p,point,current_city,borrow_city;
	double sum,cost,lost;
	for (i = 0;i < POP_SIZE;i++)
	{
		sum = 0;
		cost = 0;
        lost = 0;
		for (j = 0;j < NUM_OF_GENE;j++)
		{
			point = j / GENE_PER_POINT;
            p = point * GENE_PER_POINT + random_sort[j % GENE_PER_POINT];
			current_city = CURRENT[p % GENE_PER_POINT];
			borrow_city = BORROW[p % GENE_PER_POINT];
			//Initialize storage
			if (j % GENE_PER_POINT == 0)
				for (k = 0;k < NUM_OF_CITY;k++)
				{
					if (point > 0)
						storage[k][point] = storage[k][point - 1] + generate[k][point] - usage[k][point];
					else
						storage[k][point] = generate[k][point] - usage[k][point];
				}
			//判断是否允许借电，若电量不足则调整借电请求
			if (population[i].gene[p].borrow_amount > 0)
			{
				if (storage[borrow_city][point] >= 0)
				{
					if (storage[borrow_city][point] < population[i].gene[p].borrow_amount)
						population[i].gene[p].borrow_amount = storage[borrow_city][point];
					storage[borrow_city][point] -= population[i].gene[p].borrow_amount;
					storage[current_city][point] += population[i].gene[p].borrow_amount;
					cost += abs(population[i].gene[p].borrow_amount) * COST_PER_DIST * map[borrow_city][current_city];
				}
			}
			else if (population[i].gene[p].borrow_amount < 0)
			{
				if (storage[current_city][point] >= 0)
				{
					if (storage[current_city][point] < -population[i].gene[p].borrow_amount)
						population[i].gene[p].borrow_amount = -storage[current_city][point];
					storage[borrow_city][point] -= population[i].gene[p].borrow_amount;
					storage[current_city][point] += population[i].gene[p].borrow_amount;
					cost += abs(population[i].gene[p].borrow_amount) * COST_PER_DIST * map[borrow_city][current_city];
				}
			}
			//维护实际用电量并计算适应值
			if ((j + 1) % GENE_PER_POINT == 0)
				for (k = 0;k < NUM_OF_CITY;k++)
				{
					if (storage[k][point] < 0)
					{
						actual[k][point] = usage[k][point] + storage[k][point];
                        lost += storage[k][point] * LOST_PER_LACK;
						storage[k][point] = 0;
					}
					else actual[k][point] = usage[k][point];
					sum += PRICE_PER_UNIT * actual[k][point] - COST_PER_UNIT * generate[k][point];
				}
		}
		sum -= cost;
        if (lost > 0)
        {
            population[i].fitness = lost + 100;
            population[i].punish = sum + lost;
        }
        else
        {
            population[i].fitness = sum + lost;
            population[i].punish = lost;
        }
	}
}

void keep_best()
{
	int mem,i;
	best_record = 0;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		if (population[mem].fitness > population[POP_SIZE].fitness)
		{
			best_record = mem;
			population[POP_SIZE].fitness = population[mem].fitness;
		}
	}
    population[POP_SIZE] = population[best_record];
}

void elitist()
{
	int i;
	double best,worst;
	best = population[0].fitness;
	worst = population[0].fitness;
	for (i = 0;i < POP_SIZE - 1;i++)
	{
		if (population[i].fitness > population[i + 1].fitness)
		{
			if (population[i].fitness >= best)
			{
				best = population[i].fitness;
				best_mem = i;
			}
			if (population[i + 1].fitness <= worst)
			{
				worst = population[i + 1].fitness;
				worst_mem = i + 1;
			}
		}
		else
		{
			if (population[i].fitness <= worst)
			{
				worst = population[i].fitness;
				worst_mem = i;
			}
			if (population[i + 1].fitness >= best)
			{
				best = population[i + 1].fitness;
				best_mem = i + 1;
			}
		}
	}
	if (best > population[POP_SIZE].fitness)
	{
		for (i = 0;i < NUM_OF_GENE;i++)
			population[POP_SIZE].gene[i] = population[best_mem].gene[i];
		population[POP_SIZE].fitness = population[best_mem].fitness;
        population[POP_SIZE].punish = population[best_mem].punish;
		stable = 0;
	}
	else
	{
		for (i = 0;i < NUM_OF_GENE;i++)
			population[worst_mem].gene[i] = population[POP_SIZE].gene[i];
		population[worst_mem].fitness = population[POP_SIZE].fitness;
        population[worst_mem].punish = population[POP_SIZE].punish;
		stable++;
	}
}

void select()
{
	int first = 0,mem,i,j,one;
	double sum = 0,p,x;
	for (mem = 0;mem < POP_SIZE;mem++)
		sum += population[mem].fitness;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		population[mem].rfitness = population[mem].fitness / sum;
		population[mem].chosen = false;
	}
	population[0].cfitness = population[0].rfitness;
	for (mem = 1;mem < POP_SIZE;mem++)
		population[mem].cfitness = population[mem - 1].cfitness + population[mem].rfitness;
	/*for (i = 0;i < POP_SIZE;i++)
	{
		p = rand() % 1000 / 1000.0;
		if (p < population[0].cfitness)
			newpopulation[i] = population[0];
		else
			for (j = 0;j < POP_SIZE;j++)
				if (p >= population[j].cfitness && p < population[j + 1].cfitness)
					newpopulation[i] = population[j + 1];
	}*/
	for (i = 0;i < MODGA_R;i++)
	{
		p = rand() % 1000 / 1000.0;
		if (p < population[0].cfitness)
			newpopulation[i] = population[0];
		else
			for (j = 0;j < POP_SIZE;j++)
				if (p >= population[j].cfitness && p < population[j + 1].cfitness)
					newpopulation[i] = population[j + 1];
	}
	for (mem = 0;mem < MODGA_R;mem++)
	{
		x = rand() % 1000 / 1000.0;
		if (x < P_CROSSOVER)
		{
			first++;
			if (first % 2 == 0)
				newpopulation[mem] = newpopulation[compete(one,mem)];
			else one = mem;
		}
	}
	for (i = MODGA_R;i < POP_SIZE;i++)
	{
		p = rand() % 1000 / 1000.0;
		if (p < population[0].cfitness)
		{
			if (!population[0].chosen)
			{
				newpopulation[i] = population[0];
				population[0].chosen = true;
			}
			else i--;
		}
		else
			for (j = 0;j < POP_SIZE;j++)
				if (p >= population[j].cfitness && p < population[j + 1].cfitness)
				{
					if (!population[j + 1].chosen)
					{
						newpopulation[i] = population[j + 1];
						population[j + 1].chosen = true;
					}
					else i--;
				}
	}
	for (i = 0;i < POP_SIZE;i++)
	{
		population[i] = newpopulation[i];
		swirl(i,newpopulation[i]);
	}
}

void crossover()
{
	int first = 0,mem,one;
	double x;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		x = rand() % 1000 / 1000.0;
		if (x < P_CROSSOVER)
		{
			first++;
			if (first % 2 == 0)
				crossover(one,mem);
			else one = mem;
		}
	}
}

void crossover(int one,int two)
{
	int i,a,b;
	a = rand() % (NUM_OF_GENE - 1) + 1;
	b = rand() % (NUM_OF_GENE - a) + a;
	for (i = a;i < b;i++)
	{
		swap(population[one].gene[i],population[two].gene[i]);
		swap(population[one].lower[i],population[two].lower[i]);
		swap(population[one].upper[i],population[two].upper[i]);
	}
}

int compete(int one,int two)
{
	double p;
	if (population[one].fitness < population[two].fitness)
		return one;
	else
	{
		p = rand() % 1000 / 1000.0;
		if (p > 1 / (1 + exp((population[one].fitness + population[two].fitness) / generation)))
			return two;
		else return one;
	}
}

void mutate()
{
	int i,j;
	double x,lbound,ubound;
	for (i = 0;i < POP_SIZE;i++)
		for (j = 0;j < NUM_OF_GENE;j++)
		{
			x = (rand() % 1000) / 1000.0;
			if (x < P_MUTATION * adjust[generation])
			{
				//mutate borrow_amount
				lbound = population[i].lower[j];
				ubound = population[i].upper[j];
				population[i].gene[j].borrow_amount = random_value(lbound,ubound);
			}
		}
}

void rollback()
{
	int i,j;
	for (i = 0;i < MODGA_R / 2;i++)
	{
		j = rand() % POP_SIZE;
		swirl(j,population[POP_SIZE]);
	}
	for (i = POP_SIZE - MODGA_R / 2;i < POP_SIZE;i++)
	{
		j = rand() % POP_SIZE;
		swirl(j,history_best);
	}
    for (i = POP_SIZE / 2;i < POP_SIZE;i++)
    {
        j = rand() % POP_SIZE;
        substitute(j,POP_SIZE);
    }
	for (i = 0;i < POP_SIZE / 2;i++)
	{
		j = rand() % POP_SIZE;
		crossover(j,POP_SIZE);
	}
	population[POP_SIZE] = population[0];
}

void swirl(int index,Genetype & base)
{
	int i;
	double p;
	for (i = 0;i < NUM_OF_GENE;i++)
	{
		p = random_value(-ROLL_BOUND,ROLL_BOUND);
		population[index].gene[i].borrow_amount = base.gene[i].borrow_amount + p;
		if (population[index].gene[i].borrow_amount > population[index].upper[i])
			population[index].gene[i].borrow_amount = population[index].upper[i];
		if (population[index].gene[i].borrow_amount < population[index].lower[i])
			population[index].gene[i].borrow_amount = population[index].lower[i];
	}
}

void substitute()
{
	int first = 0,mem,one;
	double x;
	for (mem = 0;mem < POP_SIZE;mem++)
	{
		x = rand() % 1000 / 1000.0;
		if (x < P_CROSSOVER)
		{
			first++;
			if (first % 2 == 0)
				substitute(one,mem);
			else one = mem;
		}
	}
}

void substitute(int one,int two)
{
	int i,a,b;
    double x;
	for (i = 0;i < NUM_OF_GENE;i++)
	{
        x = rand() % 1000 / 1000.0;
        if (x < P_SUBSTITUTE)
        {
            population[one].gene[i] = population[two].gene[i];
            population[one].lower[i] = population[two].lower[i];
            population[one].upper[i] = population[two].upper[i];
        }
	}
}

void report(int cnt)
{
	ofstream fout(OUTPUT_FILE_NAME + to_string(cnt) + ".txt",ios::app);
	/*fout << generation << "\n";
	for (int i = 0;i < NUM_OF_CITY;i++)
	{
		fout << i + 1 << ": \n";
		for (int j = i;j < NUM_OF_GENE;j += NUM_OF_CITY)
			fout << population[POP_SIZE].gene[j].borrow_city + 1 << " " << population[POP_SIZE].gene[j].borrow_amount << "\n";
	}*/
    if (population[POP_SIZE].fitness >= population[POP_SIZE].punish)
        fout << population[POP_SIZE].fitness - population[POP_SIZE].punish << " " << population[POP_SIZE].punish << "\n";
    else fout << -population[POP_SIZE].fitness + population[POP_SIZE].punish + 100 << " " << population[POP_SIZE].fitness - 100 << "\n";
	fout.close();
}

void trace(int cnt)
{
	ofstream fout(OUTPUT_FILE_NAME + to_string(cnt) + ".txt",ios::app);
	/*fout << generation << "\n";
	for (int i = 0;i < NUM_OF_CITY;i++)
	{
		fout << i + 1 << ": \n";
		for (int j = i;j < NUM_OF_GENE;j += NUM_OF_CITY)
			fout << population[POP_SIZE].gene[j].borrow_city + 1 << " " << population[POP_SIZE].gene[j].borrow_amount << "\n";
	}
	for (int i = 0;i < NUM_OF_GENE;i++)
		fout << population[POP_SIZE].gene[i].borrow_amount << "\n";*/
    if (population[POP_SIZE].fitness > population[POP_SIZE].punish)
        fout << population[POP_SIZE].fitness - population[POP_SIZE].punish << " " << population[POP_SIZE].punish << " " << history_best.fitness - history_best.punish << " " << history_best.punish << "\n";
    else fout << -population[POP_SIZE].fitness + population[POP_SIZE].punish + 100 << " " << population[POP_SIZE].fitness - 100 << " " << -history_best.fitness + history_best.punish + 100 << " " << history_best.fitness - 100<< "\n";
	fout.close();
}

int main()
{
	int count = 0;
	string progress_bar = "";
	double complete_percent = 0;
    double border = 0;
	for (int i = 1;i <= MAX_GENERATION;i++)
		adjust[i] = 1 + exp(float(-i));
	srand((unsigned int)time(0));
    if (!VERSION)
      system("CLS");
    else printf("\033[2J\033[0;0H");
    printf("0.00%%\n");
	//init_file();
	while (count < MAX_LOOP_TIME)
	{
        border = (double)(count + 1) / (double)MAX_LOOP_TIME * 100;
		generation = 0;
		initialize(1);
		evaluate();
		keep_best();
		while (generation++ < MAX_GENERATION)
		{
			elitist();
			select();
			report(count);
            if (VERSION)
            {
                printf("\033[1A\033[K%s%.2f%%\n",progress_bar.c_str(),complete_percent + (double)generation / (double)MAX_GENERATION * 100 / (double)MAX_LOOP_TIME);
            }
            crossover();
			mutate();
			if (stable >= MAX_STABLE)
				rollback();
			if (population[POP_SIZE].fitness > history_best.fitness)
				history_best = population[POP_SIZE];
			evaluate();
		}
		trace(MAX_LOOP_TIME);
		progress_bar.append(">");
		complete_percent = border;
        if (!VERSION)
        {
            system("CLS");
    		cout << progress_bar << fixed << setprecision(2) << complete_percent << "%\n";
        }
        else
        {
            printf("\033[1A\033[K%s%.2f%%\n",progress_bar.c_str(),complete_percent);
        }
		count++;
	}
	ofstream ftotal("ga_total.txt");
	int current_city,timepoint,j;
	for (j = 0;j < NUM_OF_GENE;j++)
	{
		current_city = j % NUM_OF_CITY + 1;
		timepoint = j / NUM_OF_CITY;
		ftotal << history_best.gene[j].borrow_amount << endl;
	}
	cout << "Completed.\n";
	return 0;
}
