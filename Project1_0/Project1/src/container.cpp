// This library is free and distributed under
// Mozilla Public License Version 2.0.

// This example is demonstrating how to use 
// genetic operations which belong to a class object.


#include <string>
#include <iostream>
#include <fstream>
#include <functional>
#include "openGA.hpp"
#include <sstream>
#include <iomanip>
#include <map>
#include "math3d.h"
#include "box.h"



//struct MySolution
//{
//	double x;
//	double y;
//
//	std::string to_string() const
//	{
//		return 
//			"{x:"+std::to_string(x)+
//			", y:"+std::to_string(y)+
//			"}";
//	}
//};

using namespace Claudette::Internal;

struct MySolution

{
	std::vector<int> boxIndexes;
	std::string to_string() const
	{
		std::ostringstream out;
		out << "{";
		for (unsigned long i = 0; i < boxIndexes.size(); i++)
			out << (i ? "," : "") << std::setprecision(5) << boxIndexes[i];
		out << "}";
		return out.str();
	}
};

struct MyMiddleCost
{
	// This is where the results of simulation
	// is stored but not yet finalized.	
	double resultVolume = 0.0;
	double volumeRatio;

};

typedef EA::Genetic<MySolution,MyMiddleCost> GA_Type;
typedef EA::GenerationType<MySolution,MyMiddleCost> Generation_Type;
class Dimensions
{
public:
	int w, h, d;
	int volume;
	 Dimensions()
	{
		this->volume = 0;
	}
	 Dimensions(int w, int h,int d)
	{
		this->w = w;
		this->h = h;
		this->d = d;
		this->volume = this->w * this->d * this->h;
	}
	void calculateVolume()
	{
		this->volume = this->w * this->d * this->h;
	}

};

class BoxCapsule
{
public:
	Dimensions dims;
	int	x = 0, y = 0, z = 0;
	Claudette::Internal::Box coreBox;

	BoxCapsule()
	{

	}
	
	BoxCapsule(Dimensions d)
	{
		this->dims = d;
		Claudette::Internal::Box coreBox((float)-999990.0, (float)-999990.0, (float)-99990.0, (float)d.w, (float)d.d,(float)d.h )  ;
		this->coreBox = coreBox;
	/*	this->dims.w = d.w;
		this->dims.h = d.h;
		this->dims.d = d.d;*/
	}

};


class Container
{
public:
	std::ofstream output_file;
	int countOfBox;
	std::vector<int> coreIndexList;
	std::vector<int> randomIndexList;
	Dimensions cDims;
	std::map<int, BoxCapsule> BoxList;
	Claudette::Internal::Box containerZ;
	Claudette::Internal::Box containerX;
	Claudette::Internal::Box containerY;

	double bestTotalCost;
	int allBestGeneration_number;
	MySolution allBestGenes;

	Container(std::string output_location, Dimensions d, std::map<int, BoxCapsule> BoxList)
	{
		this->countOfBox = BoxList.size();
		this->BoxList = BoxList;
		this->cDims = d;
		Claudette::Internal::Box containerZ((float)-10, (float)-10, (float)d.h+1, (float)d.w+20, (float)d.d+20, (float)d.h+20);
		Claudette::Internal::Box containerX((float)d.w+1, (float)0, (float)0, (float)d.w + 20, (float)d.d + 20, (float)d.h + 20);
		Claudette::Internal::Box containerY((float)0, (float)d.d+1, (float)0, (float)d.w + 20, (float)d.d + 20, (float)d.h + 20);
		this->containerZ = containerZ;
		this->containerX = containerX;
		this->containerY = containerY;

		for (int i = 0; i < this->countOfBox; i++)
		{
			coreIndexList.push_back(i);
		}
		output_file.open(output_location);
		output_file<<"step"<<"\t"<<"x_best"<<"\t"<<"y_best"<<"\t"<<"cost_avg"<<"\t"<<"cost_best"<<"\n";
	}

	~Container()
	{
		output_file.close();
	}

	

	void init_genes(MySolution& p,const std::function<double(void)> &rnd01)
	{
		//static int m = 1;
		//std::cout << m << " .Gen olusturuldu " << std::endl;
		//m++;
		std::copy(this->coreIndexList.begin(), this->coreIndexList.end(),
			std::back_inserter(p.boxIndexes));
		//p.boxIndexes.swap(this->coreIndexList);
		std::random_device rd;
		std::mt19937 g(rd());

		std::shuffle(p.boxIndexes.begin(), p.boxIndexes.end(), g);
		int u = p.boxIndexes[3];
		int h = p.boxIndexes.at(1);
		
	}
	bool IsBoxIntersectAny(const MySolution& p, const BoxCapsule& input)
	{
		for (unsigned long i = 0; i < p.boxIndexes.size(); i++)
		{
			BoxCapsule currentBox = this->BoxList[p.boxIndexes[i]];
			if (input.coreBox == currentBox.coreBox)
			{
				continue;
			}
			else if (currentBox.coreBox.intersect(input.coreBox))
			{
				return true;
			}
			
			
		}
		return false;
	}

	bool eval_solution(
		const MySolution& p,
		MyMiddleCost &c)
	{		
		
		int result = 0;
		int a = p.boxIndexes.at(0);
		int b = p.boxIndexes.at(1);
		Claudette::Internal::Vector3D localOrigin;
		Claudette::Internal::Vector3D generalOrigin(0.0,0.0,0.0);
		Claudette::Internal::Box localBoundingBox(0, 0, 0, 0, 0, 0);
		Claudette::Internal::Box generalBoundingBox(0, 0, 0, 0, 0, 0);
		Matrix3D m1(1, 1, 1, 1,
			1, 1, 1, 1,
			1, 1, 1, 1,
			1, 1, 1, 1) ;
		const Matrix3D m0(0,0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0 );
		RotationState rs( m0);
		
		for (unsigned long i = 0; i < p.boxIndexes.size(); i++)
		{

			BoxCapsule currentBox = this->BoxList[p.boxIndexes[i]];
			
			currentBox.coreBox.m_Pos = localBoundingBox.getPosition() + Vector3D(0,0,localBoundingBox.getSize().z);
			
			if (currentBox.coreBox.intersect(this->containerZ))
			{
				currentBox.coreBox.m_Pos = localBoundingBox.getPosition() + Vector3D(localBoundingBox.getSize().x,0, 0);
				if (currentBox.coreBox.intersect(this->containerX))
				{
					currentBox.coreBox.m_Pos = generalBoundingBox.getPosition() + Vector3D(0, generalBoundingBox.getSize().y, 0);
					if (currentBox.coreBox.intersect(this->containerY))
					{
						//settleState = false;
						break;
					}
					else
					{
						//settleState = true;
						localBoundingBox = currentBox.coreBox;
					}
				}
				else
				{
					localBoundingBox.m_Size = currentBox.coreBox.m_Size;					
					localBoundingBox.m_Pos = currentBox.coreBox.m_Pos;
					generalBoundingBox.m_Size = currentBox.coreBox.m_Pos + currentBox.coreBox.m_Size;
				}
			}
			else
			{

				if ((localBoundingBox.getSize().y > 0 && (localBoundingBox.getSize().y * 2 < currentBox.coreBox.getSize().y)) ||
					(localBoundingBox.getSize().x > 0 && (localBoundingBox.getSize().x * 2 < currentBox.coreBox.getSize().x)))
				{
					//break;
				}
					localBoundingBox.m_Size = Vector3D(localBoundingBox.getSize().x > currentBox.coreBox.getSize().x ? localBoundingBox.getSize().x : currentBox.coreBox.getSize().x,
					localBoundingBox.getSize().y > currentBox.coreBox.getSize().y ? localBoundingBox.getSize().y : currentBox.coreBox.getSize().y,
					localBoundingBox.getSize().z +currentBox.coreBox.getSize().z);
					generalBoundingBox.m_Size = currentBox.coreBox.m_Pos + currentBox.coreBox.m_Size;

				 
			}

			if (!IsBoxIntersectAny(p,currentBox))
			{
				lastPlaced = i;
				//localBoundingBox.m_Size += Vector3D(currentBox.dims.w, currentBox.dims.d, currentBox.dims.h);
				//generalBoundingBox.m_Size += Vector3D(currentBox.dims.w, currentBox.dims.d, currentBox.dims.h);
				result += this->BoxList[p.boxIndexes[i]].dims.volume;
			}
			
			
		}
		c.resultVolume = result;
		c.volumeRatio = ( (double) result / (double)this->cDims.volume) ;
		return true;

	}

	MySolution mutateCore(
		const MySolution& X_base,
		const std::function<double(void)> &rnd01
		)
	{
		MySolution X_new;
		int randomIndex0 = this->countOfBox * rnd01();
		int randomIndex1 = this->countOfBox * rnd01();
		for (unsigned long i = 0; i < X_base.boxIndexes.size(); i++)
		{
			X_new.boxIndexes.push_back(X_base.boxIndexes[i]);
		}
		X_new.boxIndexes[randomIndex0] = X_base.boxIndexes.at(randomIndex1);
		X_new.boxIndexes[randomIndex1] = X_base.boxIndexes.at(randomIndex0);
		return X_new;
	}

	MySolution mutate(
		const MySolution& X_base,
		const std::function<double(void)>& rnd01,
		double shrink_scale)
	{
		MySolution X_new;
		X_new = X_base;
		for (unsigned long i = 0; i < this->BoxList.size()* shrink_scale; i++)
		{
			X_new = mutateCore(X_new, rnd01);
		}


		return X_new;
	}

	MySolution crossover(
		const MySolution& X1,
		const MySolution& X2,
		const std::function<double(void)> &rnd01)
	{
		MySolution  X_new;
		int randomIndex = this->countOfBox * rnd01();
		int randomCount = ((this->countOfBox- randomIndex) * rnd01()) * 0.5;
		std::vector<int> dst = X1.boxIndexes;
		std::vector<int> dstNew;
		std::vector<int>::iterator it;	
		
		std::copy(dst.begin()+ randomIndex, dst.begin() + randomIndex+ randomCount, std::back_inserter(dstNew));
			for(unsigned long i=0;i<X2.boxIndexes.size() ;i++)
			{			
				it = std::find(dstNew.begin(), dstNew.end(), X2.boxIndexes[i]);
				if (it == dstNew.end())
				{
					X_new.boxIndexes.push_back(X2.boxIndexes[i]);
				}			

			}
			X_new.boxIndexes.insert(X_new.boxIndexes.begin() + randomIndex, dstNew.begin(), dstNew.end());	
		
		return X_new;
	}

	double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X)
	{
		// finalize the cost
		//double cost1,cost2;
		//cost1=X.middle_costs.cost_distance2;
		//cost2=X.middle_costs.cost_sqsin;
		return X.middle_costs.volumeRatio ;
	}

	void SO_report_generation(		
		int generation_number,
		const EA::GenerationType<MySolution,MyMiddleCost> &last_generation,
		const MySolution& best_genes)
	{



		if (last_generation.best_total_cost  > this->bestTotalCost)
		{
			this->bestTotalCost = last_generation.best_total_cost;
			this->allBestGeneration_number = generation_number;
			this->allBestGenes= best_genes;
		}

		std::cout << " Nesil [" << generation_number << "] " << std::endl;
		std::cout		
			<<" En iyi="<<last_generation.best_total_cost <<", "		
			<<" En iyi gen=("<<best_genes.to_string()<<")"<<", "
			<<" islem suresi ="<<last_generation.exe_time
			<<std::endl;

		std::cout << " " << std::endl;

		output_file
			<<generation_number<<"\t"
			<<best_genes.to_string()<<"\t"
			<<best_genes.to_string() <<"\t"
			<<last_generation.average_cost<<"\t"
			<<last_generation.best_total_cost<<"\n";
	}
};

std::map<int, BoxCapsule> initialBoxlistManual()
{
	std::map<int, BoxCapsule> boxList = {
	   {0, BoxCapsule(Dimensions(20, 21, 22))},
	   {1, BoxCapsule(Dimensions(20, 21, 22))},
	   {2, BoxCapsule(Dimensions(40, 21, 22))},
	   {3, BoxCapsule(Dimensions(50, 21, 22))},
	   {4, BoxCapsule(Dimensions(60, 21, 22))},
	   {5, BoxCapsule(Dimensions(20, 21, 22))},
	   {6, BoxCapsule(Dimensions(20, 21, 22))},
	   {7, BoxCapsule(Dimensions(40, 21, 22))},
	   {8, BoxCapsule(Dimensions(50, 21, 22))},
	   {9, BoxCapsule(Dimensions(60, 21, 22))},
	   {10, BoxCapsule(Dimensions(60, 21, 22))},
	};
	return	boxList;
}

std::map<int, BoxCapsule> initialBoxlistSimple()
{
	std::map<int, BoxCapsule> boxList;

	// https://matematikdelisi.com/Orta/Sinif8/Konu/EBOBEKOK7/img/kuplerden-dikdortgen-prizma-olusturma.png
	// basit kontrol
	for (int i = 0; i < 17; i++)
	{
		boxList.insert(std::pair<int, BoxCapsule>(i, BoxCapsule(Dimensions(10, 10, 10))));
	}
	//for (int i = 0; i < 1; i++)
	//{
	//	boxList.insert(std::pair<int, BoxCapsule>(15, BoxCapsule(Dimensions(10, 10, 10))));
	//}

	

	
	 return boxList;

}

std::map<int, BoxCapsule> initialTHPackBoxlist()
{
	std::map<int, BoxCapsule> boxList;


	for (int i = 0; i < 39; i++)
	{
		boxList.insert(std::pair<int, BoxCapsule>(i, BoxCapsule(Dimensions(108, 76, 30))));
	}


	for (int i = 39; i < (39+32); i++)
	{
		boxList.insert(std::pair<int, BoxCapsule>(i, BoxCapsule(Dimensions(110, 43, 25))));
	}

	for (int i = (39 + 32); i < (39 + 32+38); i++)
	{
		boxList.insert(std::pair<int, BoxCapsule>(i, BoxCapsule(Dimensions(92, 81, 55))));
	}


	return boxList;

}

int main()
{
	char n;
	
	do
	{
		std::cout << '\n' << "Baslamak icin bir tusa basin..";
	} while (std::cin.get() != '\n');
	
	std::map<int, BoxCapsule> boxList = initialTHPackBoxlist();
	//boxList[1] = Box(Dimensions(20, 20, 20));


	// basit kontrol 
	//Container myobject("./bin/result.txt",Dimensions(40,20,20), boxList);
	
	//thpack
	Container myobject("./bin/result.txt", Dimensions(587, 233, 220), boxList);

	EA::Chronometer timer;
	timer.tic();

	using std::bind;
	using std::placeholders::_1;
	using std::placeholders::_2;
	using std::placeholders::_3;


	GA_Type ga_obj;
	ga_obj.problem_mode= EA::GA_MODE::SOGA;
	ga_obj.multi_threading=true;
	ga_obj.idle_delay_us=1; 
	ga_obj.verbose=false;
	ga_obj.population=500;
	ga_obj.generation_max=100;
	ga_obj.calculate_SO_total_fitness=bind(&Container::calculate_SO_total_fitness, &myobject, _1);
	ga_obj.init_genes= bind(&Container::init_genes, &myobject, _1, _2);
	ga_obj.eval_solution= bind(&Container::eval_solution, &myobject, _1, _2);
	ga_obj.mutate= bind(&Container::mutate, &myobject, _1, _2, _3);
	ga_obj.crossover= bind(&Container::crossover, &myobject, _1, _2, _3);
	ga_obj.SO_report_generation= bind(&Container::SO_report_generation, &myobject, _1, _2, _3);
	ga_obj.average_stall_max = 5000;
	ga_obj.best_stall_max=5000;
	ga_obj.elite_count=10;
	ga_obj.crossover_fraction=0.5;
	ga_obj.mutation_rate=0.1;
	ga_obj.solve();


	std::cout <<  "---EN IYI SONUC---" << std::endl;
	std::cout << "  Nesil [" << myobject.allBestGeneration_number << "] " << std::endl;
	std::cout
		<< " En iyi=" << myobject.bestTotalCost << ", "
		<< " En iyi gen=(" << myobject.allBestGenes.to_string() << ")" << ", "

		<< std::endl;

	std::cout << " " << std::endl;

	std::cout<<"Problem optimizasyonu suresi: "<<timer.toc()<<" saniye."<<std::endl;
	

	return 0;
}

