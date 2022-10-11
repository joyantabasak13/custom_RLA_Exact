/*
 * RLA_CL.cpp
 *
 *  Created on: Apr 4, 2014
 *      Author: abdullah-al-mamun
 */


// library
#include <iostream>
#include <fstream>
#include <libxml/xpath.h>
#include <libxml/parser.h>
#include </opt/homebrew/Cellar/boost/1.79.0_2/include/boost/algorithm/string.hpp>
#include <algorithm>
#include <iterator>
#include <vector>
#include </opt/homebrew/Cellar/boost/1.79.0_2/include/boost/lexical_cast.hpp>
#include <ctype.h>
#include <math.h>
#include <cmath>
#include <utility>

#include "Cluster.h"
#include "Common.h"
#include "SortMSD.h"

// namespace
using namespace std;
using namespace boost;

// MACRO
#define RECORD_TOTAL		500
#define THRESHOLD_MULTI		100
#define NUMBER_LMER			3
#define ROOT_ID				0
#define LINE_LENGTH_MIN		10
#define DATA_INVALID		-1
#define TYPE_EDIT			0
#define TYPE_REVERSAL		1
#define TYPE_TRUNC			2
#define TYPE_TOTAL			3

#define OUTPUT_FINAL		0
#define OUTPUT_SINGLE		1

// structure
typedef struct
{
	int ind, clusterInd;
}RecordPacket;



vector<int> getWeight();
vector<vector<int> > getIndexDataSet();
vector<string> getInputFileNameList();
int getInputThreshold();
vector<vector<int> > getBlockField();
vector<int> getPriorityField();
vector<vector<vector<int> > > getInputComparisonPara();
string getOutputDir();

void readDataFromFile(vector<string>& fileNameArr);

void clusterData();
void sortData();
void findExactClusterSL();
void findExactClusterEach();
void findExactCluster();
void appendExactInd(int pos, vector<int>& tempArr);
void findExactClusterPart(int i, int j);
void findApproxCluster(vector<int>& rootArr);
void findFinalCluster(vector<int>& rootArr);

void createBlock(int indBlockField, int lmerUsed, int type, vector<vector<int> >& blockArr);
void createClusterEdgeList(vector<vector<int> >& blockArr);
void generateEdgilist(vector<int>& blockRowArr);
void findConnComp(vector<int>& parentArr);

bool isLinkageOk(vector<string>& a, vector<string>& b);
int linkage(vector<string>& a, vector<string>& b);
int compareAttLen(vector<string>& a, vector<string>& b);
int calculateDistAll(vector<string>& a, vector<string>& b);
int calculateEditDist(vector<string>& a, vector<string>& b, vector<vector<int> >& compareAtt, int threshRem);
int calculateRevDist(vector<string>& a, vector<string>& b, vector<vector<int> >& compareAtt, int threshRem);
int calculateTruncDist(vector<string>& a, vector<string>& b, vector<vector<int> >& compareAtt, int threshRem);
int calculateBasicED(string& str1, string& str2, int threshRem);
int calculateBasicED2(string& str1, string& str2, int threshRem);
string convertToLower(string& str);
int findRoot(int pointID, vector<int>& parentArr);
void makeEquivalent(int rootU, int rootV, vector<int>& parentArr, vector<int>& weightArr);

void clusterGrp(vector<int>& recordIndSingleGrpArr, vector<Cluster>& clusterSingleArr);
void generateMatrix(vector<vector<string> >& recordGrpArr, vector<vector<int> >& matArr);
void generateVector(vector<vector<int> >& matArr, vector<vector<int> >& vecArr);
void updateMatVec(vector<vector<int> > &matArr, vector<vector<int> > &vecArr, int indCluster1, int indCluster2);
void postprocessCluster(vector<Cluster>& clusterSingleArr);
bool isInCluster(Cluster& cluster, vector<string>& record);
void processUsedAttr(vector<vector<int> >& usedAttrArr, vector<int>& isInClusterArr, int indDataset, vector<vector<int> >& usedThisAttrArr);

void output(int outputType);

static const int ALPHABET_SIZE_LIST[] = {26, 10, 36, 256};

int testVal;
clock_t startT;
double readT, clusterExactT, clusterApproxT, clusterFinalT, clusterCLT, totalT, outT;
int recordTotal, recordTotalCMD;
string configFileStr, outFileAddrStr, inFileAddrStr, strSample, outDir;
int threshold;

//Ruru
int rand_theshold = 30;
bool is_rand = true;

vector<Cluster> clusterArr;
vector<string> fileNameArr;
vector<vector<string> > recordArr;
vector<int> weightArr, priorityFieldArr, recordStartIndArr;
vector<vector<int> > edgeArr, clusterExactIndArr, indexDatasetArr, clusterIndArr, blockFieldArr, clusterExactPairArr;
vector<vector<vector<int> > > attrArr;

int lenMax;
vector<pair<int, string>> strDataArr;

void doRLA(string fileAddr, int recordTotal)
{
	cout << "initRLA()" << endl;

	// read data from configuration file
	recordTotalCMD	= recordTotal;
	configFileStr	= fileAddr;

	weightArr		= getWeight();
	indexDatasetArr	= getIndexDataSet();
	fileNameArr		= getInputFileNameList();
	blockFieldArr	= getBlockField();
	priorityFieldArr= getPriorityField();
	threshold		= getInputThreshold();
	
	
	attrArr			= getInputComparisonPara();
	outDir			= getOutputDir();

	// read data from input file
	startT		= clock();
	readDataFromFile(fileNameArr);
	readT		= (double)(clock() - startT) / CLOCKS_PER_SEC;;

	cout << recordArr.size() << " records; Threshold " << threshold << endl;

	// cluster data
	clusterData();
}

// read weight parameters from config xml
vector<int> getWeight()
{
	cout << "getWeight()" << endl;

	xmlDoc *doc 				= xmlReadFile(configFileStr.c_str(), NULL, 0);
	xmlXPathContext *pathCtx 	= xmlXPathNewContext(doc);
	xmlXPathObject *pathObj 	= xmlXPathEvalExpression((xmlChar *)"/rla-config/weights", pathCtx);
	xmlXPathFreeContext(pathCtx);

	vector<int> weightArr;
	for (int i = 0; i < pathObj->nodesetval->nodeNr; ++i)
	{
		xmlNode *nodeCurr 	= pathObj->nodesetval->nodeTab[i];
		xmlNode *nodeChild 	= nodeCurr->children;

		while(nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
			{
				string weight = (char *)xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1);
				vector<string> valStrArr;
				split(valStrArr, weight, is_any_of(" "));
				for(vector<string>::iterator it = valStrArr.begin(); it != valStrArr.end(); ++it)
					weightArr.push_back(atoi((*it).c_str()));
				break;
			}
			nodeChild = nodeChild->next;
		}
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return weightArr;
}

// read index of attributes of data sets from config xml
vector<vector<int> > getIndexDataSet()
{
	cout << "getIndexDataSet" << endl;

	xmlDoc *doc 				= xmlReadFile(configFileStr.c_str(), NULL, 0);
	xmlXPathContext *pathCtx 	= xmlXPathNewContext(doc);
	xmlXPathObject *pathObj 	= xmlXPathEvalExpression((xmlChar *)"/rla-config/dataset/dataset_index", pathCtx);
	xmlXPathFreeContext(pathCtx);

	vector<vector<int> > indexArr;
	for (int i = 0; i < pathObj->nodesetval->nodeNr; ++i)
	{
		xmlNode *nodeCurr = pathObj->nodesetval->nodeTab[i];
		xmlNode *nodeChild = nodeCurr->children;

		while(nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
			{
				string weight = (char *)xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1);
				vector<string> valStrArr;
				split(valStrArr, weight, is_any_of(" "));
				vector<int> indexRowArr;
				for(vector<string>::iterator it = valStrArr.begin(); it != valStrArr.end(); ++it)
					indexRowArr.push_back(atoi((*it).c_str()));
				indexArr.push_back(indexRowArr);
				break;
			}
			nodeChild = nodeChild->next;
		}
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return indexArr;
}

// read input file names from config xml
vector<string> getInputFileNameList()
{
	cout << "getInputFileNameList" << endl;

	xmlDoc *doc					= xmlReadFile(configFileStr.c_str(), NULL, 0);
	xmlXPathContext *xpathCtx	= xmlXPathNewContext(doc);
	xmlXPathObject *xpathObj	= xmlXPathEvalExpression((xmlChar *)"/rla-config/dataset", xpathCtx);
	xmlXPathFreeContext(xpathCtx);

	vector<string> fileListArr;
	for (int i = 0; i < xpathObj->nodesetval->nodeNr; ++i)
	{
		xmlNode *nodeCurr	= xpathObj->nodesetval->nodeTab[i];
		xmlNode *nodeChild	= nodeCurr->children;
		while (nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
			{
				fileListArr.push_back((char *)xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1));
				break;
			}
			nodeChild	= nodeChild->next;
		}

	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return fileListArr;
}

// read threshold value from config xml
int getInputThreshold()
{
	cout << "getInputThreshold" << endl;

	xmlDoc *doc					= xmlReadFile(configFileStr.c_str(), NULL, 0);
	xmlXPathContext *xpathCtx	= xmlXPathNewContext(doc);
	xmlXPathObject *xpathObj	= xmlXPathEvalExpression((xmlChar *)"/rla-config/version-config-param/threshold", xpathCtx);
	xmlXPathFreeContext(xpathCtx);
cout << "getInputThreshold:" << endl;
	xmlChar *attrVal	= (xmlChar *)"";
	xmlNode *nodeCurr	= xpathObj->nodesetval->nodeTab[0];
	xmlNode *nodeChild	= nodeCurr->children;
	cout << "getInputThreshold::" << endl;

	while (nodeChild)
	{
		if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
		{
			attrVal	= xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1);
			break;
		}
		nodeChild	= nodeChild->next;
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();
cout << "getInputThreshold" << endl;
	return atoi((char *)attrVal);
}

// read block field index, length & type from config xml
vector<vector<int> > getBlockField()
{
	cout << "getBlockField" << endl;

	xmlDoc *doc 				= xmlReadFile(configFileStr.c_str(), NULL, 0);
	xmlXPathContext *pathCtx 	= xmlXPathNewContext(doc);
	xmlXPathObject *pathObj 	= xmlXPathEvalExpression((xmlChar *)"/rla-config/version-config-param/block", pathCtx);
	xmlXPathFreeContext(pathCtx);

	vector<vector<int> > blockFieldArr(3);
	for (int i = 0; i < pathObj->nodesetval->nodeNr; ++i)
	{
		xmlNode *nodeCurr = pathObj->nodesetval->nodeTab[i];
		xmlNode *nodeChild = nodeCurr->children;

		while(nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"index"))
			{
				xmlNode *valIndex	= nodeChild->children;
				while(valIndex)
				{
					if(!xmlStrcmp(valIndex->name, (xmlChar *)"value"))
					{
						string index = (char *)xmlNodeListGetString(nodeCurr->doc, valIndex->children, 0);
						vector<string> valStrArr;
						split(valStrArr, index, is_any_of(","));
						vector<int> indexRowArr;
						for(vector<string>::iterator it = valStrArr.begin(); it != valStrArr.end(); ++it)
							indexRowArr.push_back(atoi((*it).c_str()));
						blockFieldArr.at(0) = indexRowArr;
						break;
					}
					valIndex = valIndex->next;
				}
			}
			else if(!xmlStrcmp(nodeChild->name, (xmlChar *)"length"))
			{
				xmlNode *valLength	= nodeChild->children;
				while(valLength)
				{
					if(!xmlStrcmp(valLength->name, (xmlChar *)"value"))
					{
						string length = (char *)xmlNodeListGetString(nodeCurr->doc, valLength->children, 0);
						vector<string> valStrArr;
						split(valStrArr, length, is_any_of(","));
						vector<int> lengthArr;
						for(vector<string>::iterator it = valStrArr.begin(); it != valStrArr.end(); ++it)
							lengthArr.push_back(atoi((*it).c_str()));
						blockFieldArr.at(1) = lengthArr;
						break;
					}
					valLength = valLength->next;
				}
			}
			else if(!xmlStrcmp(nodeChild->name, (xmlChar *)"type"))
			{
				xmlNode *valType	= nodeChild->children;
				while (valType)
				{
					if(!xmlStrcmp(valType->name, (xmlChar *)"value"))
					{
						string type = (char *)xmlNodeListGetString(nodeCurr->doc, valType->children, 0);
						vector<string> valStrArr;
						split(valStrArr, type, is_any_of(","));
						vector<int> typeArr;
						for(vector<string>::iterator it = valStrArr.begin(); it != valStrArr.end(); ++it)
							typeArr.push_back(atoi((*it).c_str()));
						blockFieldArr.at(2) = typeArr;
						break;
					}
					valType = valType->next;
				}

			}
			nodeChild = nodeChild->next;
		}
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return blockFieldArr;
}

// read priority field index from config xml
vector<int> getPriorityField()
{
	cout << "getPriorityField" << endl;

	xmlDoc *doc 				= xmlReadFile(configFileStr.c_str(), NULL, 0);
	xmlXPathContext *pathCtx 	= xmlXPathNewContext(doc);
	xmlXPathObject *pathObj 	= xmlXPathEvalExpression((xmlChar *)"/rla-config/version-config-param/priority", pathCtx);
	xmlXPathFreeContext(pathCtx);

	vector<int> priorityFieldArr;
	for (int i = 0; i < pathObj->nodesetval->nodeNr; ++i)
	{
		xmlNode *nodeCurr = pathObj->nodesetval->nodeTab[i];
		xmlNode *nodeChild = nodeCurr->children;

		while(nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
			{
				string priority = (char *)xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1);
				vector<string> valStrArr;
				split(valStrArr, priority, is_any_of(","));
				for(vector<string>::iterator it = valStrArr.begin(); it != valStrArr.end(); ++it)
					priorityFieldArr.push_back(atoi((*it).c_str()));
				break;
			}
			nodeChild = nodeChild->next;
		}
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return priorityFieldArr;
}


// read comparison methods from config xml
vector<vector<vector<int> > > getInputComparisonPara()
{
	cout << "getInputComparisonPara" << endl;

	xmlDoc *doc 				= xmlReadFile(configFileStr.c_str(), NULL, 0);
	xmlXPathContext *pathCtx 	= xmlXPathNewContext(doc);
	xmlXPathObject *pathObj 	= xmlXPathEvalExpression((xmlChar *)"/rla-config/version-config-param/comparison", pathCtx);
	xmlXPathFreeContext(pathCtx);

	xmlChar *attrVal;
	vector<int> attrSingleArr;
	vector<vector<int> > attrTypeArr;
	vector<vector<vector<int> > > attrArr(3);
	for (unsigned int i = 0; i < pathObj->nodesetval->nodeNr; ++i)
	{
		attrSingleArr.clear();
		cout << i << "\t" << pathObj->nodesetval->nodeNr << endl;
		xmlNode *nodeCurr	= pathObj->nodesetval->nodeTab[i];
		xmlNode *nodeChild	= nodeCurr->children;
		while(nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"dist_calc_method"))
			{
				xmlNode *nodeGChild	= nodeChild->children;
				while(nodeGChild)
				{
					if(!xmlStrcmp(nodeGChild->name, (xmlChar *)"value"))
					{
						attrVal		= xmlNodeListGetString(nodeCurr->doc, nodeGChild->children, 1);
						break;
					}
					nodeGChild	= nodeGChild->next;
				}
			}

			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"comparing_attribute_indices"))
			{
				xmlNode *nodeGChild	= nodeChild->children;
				while(nodeGChild)
				{
					if(!xmlStrcmp(nodeGChild->name, (xmlChar *)"value"))
					{

						if(atoi((char *)attrVal) == 1) // edit distance
						{
							attrSingleArr.push_back(atoi((char *)xmlNodeListGetString(nodeCurr->doc, nodeGChild->children, 1))); // indices
							attrSingleArr.push_back(0);
							attrArr[0].push_back(attrSingleArr);
						}
						else if(atoi((char *)attrVal) == 2) // reversal distance
						{
							//cout << attrVal << " aastt" << endl;
							string revStr	= string((char *)xmlNodeListGetString(nodeCurr->doc, nodeGChild->children, 1));
							revStr.erase(remove(revStr.begin(), revStr.end(), ' '), revStr.end());

							vector<string> valStrArr;
							split(valStrArr, revStr, is_any_of(","));
							for (unsigned int t = 0; t < valStrArr.size(); ++t)
								attrSingleArr.push_back(atoi(valStrArr[t].c_str()));
							attrArr[1].push_back(attrSingleArr);
						}
						else if(atoi((char *)attrVal) == 3) // truncation distance
						{
							attrSingleArr.push_back(atoi((char *)xmlNodeListGetString(nodeCurr->doc, nodeGChild->children, 1))); // indices

						}

						break;
					}
					nodeGChild	= nodeGChild->next;
				}
			}

			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"truncate_count"))
			{
				xmlNode *nodeGChild	= nodeChild->children;
				while(nodeGChild)
				{
					if(!xmlStrcmp(nodeGChild->name, (xmlChar *)"value"))
					{
						if(atoi((char *)attrVal) == 3) // truncation distance
						{
							attrSingleArr.push_back(atoi((char *)xmlNodeListGetString(nodeCurr->doc, nodeGChild->children, 1))); // truncation count
							attrArr[2].push_back(attrSingleArr);

							break;
						}

					}
					nodeGChild	= nodeGChild->next;
				}
			}
			nodeChild	= nodeChild->next;
		}
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return attrArr;
}


string getOutputDir()
{
	xmlDoc *doc;
	xmlNode *nodeCurr, *nodeChild, *nodeRoot;
	xmlChar *attrVal;
	xmlXPathContext *xpathCtx;
	xmlXPathObject *xpathObj;

	doc	= xmlReadFile(configFileStr.c_str(), NULL, 0);

	nodeRoot	= xmlDocGetRootElement(doc);
	xpathCtx	= xmlXPathNewContext(doc);
	xpathObj	= xmlXPathEvalExpression((xmlChar *)"/rla-config/version-config-param/output_function/output_filename", xpathCtx);
	xmlXPathFreeContext(xpathCtx);

	attrVal		= (xmlChar *)"";
	nodeCurr	= xpathObj->nodesetval->nodeTab[0];
	nodeChild	= nodeCurr->children;
	while(nodeChild)
	{
		if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
		{
			attrVal	= xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1);
			break;
		}
		nodeChild	= nodeChild->next;
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return (char *)attrVal;
}

// read data from input files
void readDataFromFile(vector<string>& fileNameArr)
{
	// store total attributes of each data set
	vector<int> fieldSizeArr;
	unsigned int k;
	for (unsigned int i = 0; i < fileNameArr.size(); ++i)
	{
		k = 0;
		for (unsigned int j = 0; j < indexDatasetArr.at(i).size(); ++j)
		{
			if (indexDatasetArr.at(i).at(j) < 0)
				continue;
			++k;
		}
		fieldSizeArr.push_back(k);
	}

	// read data
	unsigned int count;
	int recordPerFile = ceil(recordTotalCMD / fileNameArr.size());
	for (unsigned int i = 0; i < fileNameArr.size(); ++i)
	{
		recordStartIndArr.push_back(recordArr.size());
		count = 0;
		cout << fileNameArr[i] << endl;
		ifstream inFile;
		inFile.open(fileNameArr[i].c_str(), ifstream::in);
		string lineStr;
		vector<string> rowStrArr;
		if (inFile.good())
		{
			while (getline(inFile, lineStr))
			{
				if(lineStr.length() >= LINE_LENGTH_MIN)
				{
					rowStrArr.clear();
					vector<string> tempStrArr;
					split(tempStrArr, lineStr, is_any_of(","));
					for (unsigned int t = 0; t < tempStrArr.size(); ++t)
					{
						tempStrArr[t].erase(std::remove_if(tempStrArr[t].begin(), tempStrArr[t].end(), ::isspace), tempStrArr[t].end());
						rowStrArr.push_back(convertToLower(tempStrArr[t]));
					}

					for (unsigned int j = rowStrArr.size(); j < fieldSizeArr[i]; ++j)
						rowStrArr.push_back("");

					rowStrArr.push_back(fileNameArr[i]);
					rowStrArr.push_back(lexical_cast<string>(i));
					recordArr.push_back(rowStrArr);

					if (++count >= recordPerFile)
						break;
				}
			}
		}

		inFile.close();
	}

	recordTotal = recordArr.size();
	
	
	// concat common attributes for sorting
	string strSample(50, '0');	
	int attrCommonTotal = attrArr[TYPE_EDIT].size();

	lenMax		= 0;
	strDataArr.resize(recordTotal);
	
	for (int i = 0; i < recordTotal; ++i) {
		string s;
		for (int j = 0; j < attrCommonTotal; ++j) 
			s += recordArr[i][attrArr[TYPE_EDIT][j][0]] + "#";
		strDataArr[i]	= make_pair(i, s);
		if(s.length() > lenMax)
			lenMax		= s.length();
		//cout << i << ":" << s << " len: " << s.length() << endl;
	}
	for (int i = 0; i < recordTotal; ++i) {
		int lenDiff		= lenMax - strDataArr[i].second.length();
		if(lenDiff > 0)
			strDataArr[i].second	+= strSample.substr(0, lenDiff);
	}
}

void findExactClusterEach()
{
	for (int i = 0; i < recordTotal; ++i)
	{
		vector<int> rowArr;
		rowArr.push_back(i);
		clusterExactIndArr.push_back(rowArr);
	}
}

/*
 * sort(radix sort) records using some attributes
 */
void sortData() {
	//cout << lenMax << " sorting " << recordTotal << endl;
	vector<pair<int, string>> tempArr(recordTotal);
	
	for (int i = lenMax - 1; i >= 0; --i) {
		vector<int> countArr(256, 0);
		
		for (int j = 0; j < recordTotal; ++j) {
			assert(strDataArr[j].second.length() == lenMax);
			countArr[(strDataArr[j].second)[i]]++;
		}
		
		for (int k = 1; k < 256; ++k)
			countArr[k]	+= countArr[k - 1];
		
		for (int j = recordTotal - 1; j >= 0; --j)
			tempArr[--countArr[(strDataArr[j].second)[i]]]	= strDataArr[j];
		
		for (int j = 0; j < recordTotal; ++j)
			strDataArr[j]	= tempArr[j];
	}
}


void PrintData() {
	cout << "printing data " << endl;
	for (int i = 0; i < strDataArr.size(); ++i) {
		cout << strDataArr[i].first << "\t" << strDataArr[i].second << endl;
	}
}

/*
 * find clusters with no errors by grouping records using exact match 
 */
void findExactClusterSL() {
	sortData(); // sort records using some attributes as values (radix sort)
	
	//PrintData();
	
	vector<int> clusterRowArr;
	clusterRowArr.push_back(strDataArr[0].first);
	
	for (int i = 1; i < recordTotal; ++i) {
		if(strDataArr[i].second.compare(strDataArr[i - 1].second) == 0)
			clusterRowArr.push_back(strDataArr[i].first);
		else {
			clusterExactIndArr.push_back(clusterRowArr);
			clusterRowArr.clear();
			clusterRowArr.push_back(strDataArr[i].first);
		}
	}
	clusterExactIndArr.push_back(clusterRowArr);
	
	cout << "total exact clusters: " << clusterExactIndArr.size() << endl;
}




// cluster data
void clusterData()
{
	cout << "clusterData()" << endl;

	startT			= clock();

	// find exact clusters
	clock_t currTS1	= clock();
	//findExactClusterEach();
	//findExactCluster(); // find clusters using exact match
	findExactClusterSL(); // find clusters using exact match
	cout << "# of exact clusters: " << clusterExactIndArr.size() << endl;
	clusterExactT	= (double)(clock() - currTS1) / CLOCKS_PER_SEC;

	// find approximate clusters
	clock_t currTS2	= clock();
	vector<int> rootArr;
	findApproxCluster(rootArr); // find clusters using approximate match
	clusterApproxT	= (double)(clock() - currTS2) / CLOCKS_PER_SEC;


	// find all clusters using single linkage clustering
	clock_t currTS3	= clock();
	findFinalCluster(rootArr); // combine exact and approximate match
	clusterFinalT	= (double)(clock() - currTS3) / CLOCKS_PER_SEC;

	clock_t currTS9		= clock();
	output(OUTPUT_SINGLE);
	double totalOut1T	= (double)(clock() - currTS9) / CLOCKS_PER_SEC;


	cout << "\nFinding complete linkage " << endl;

	// find clusters using complete linkage clustering
	clock_t currTS4	= clock();

	for (unsigned int i = 0, m = clusterIndArr.size(); i < m; ++i)
	{
		if (clusterIndArr[i].size() > 1)
		{
			vector<Cluster> clusterTempArr;
			clusterGrp(clusterIndArr[i], clusterTempArr);
			for (unsigned int j = 0, n = clusterTempArr.size(); j < n; ++j)
				clusterArr.push_back(clusterTempArr[j]);
			//cout << "\nClustering Done with size " << clusterArr.size() << endl;
		}
		else if(clusterIndArr[i].size() == 1)
		{
			Cluster newC;
			newC.initCluster(0, recordArr[clusterIndArr[i][0]]);
			clusterArr.push_back(newC);
		}

	}
	clusterCLT		= (double)(clock() - currTS4) / CLOCKS_PER_SEC;

	totalT			= (double)(clock() - startT) / CLOCKS_PER_SEC;

	cout << "Final Clustering Done: " << clusterArr.size() << endl;

	clock_t currTS5	= clock();
	output(OUTPUT_FINAL); // write output to a file
	outT			= (double)(clock() - currTS5) / CLOCKS_PER_SEC;

	cout << "recordTotal: " << recordTotal << " readT: " << readT << " outT: " << outT << endl;
	cout << "totalT: " << (totalT - totalOut1T) << " exactT: " << clusterExactT << " approxT: " << clusterApproxT << " finalT: " << clusterFinalT << " completeT: " << clusterCLT << endl;

}




// find exact clusters
void findExactCluster()
{
	cout << "findExactCluster()" << endl;

	clusterExactPairArr.resize(recordTotal);
	for (unsigned int i = 0; i < recordTotal; ++i)
		clusterExactPairArr[i].push_back(i);

	// find pairs of datasets such that dataset j contains dataset i fully
	for (unsigned int i = 0, m = indexDatasetArr.size(); i < m; ++i)
	{
		for (unsigned int j = 0; j < m; ++j)
		{
			if (i == j)	continue;
			bool sameF	= true; // if both have same indices of attributes
			unsigned int k = 0, n = indexDatasetArr.at(i).size();
			for (; k < n; ++k)
			{
				if (indexDatasetArr[i][k] > -1)
				{
					if (indexDatasetArr[j][k] == -1)	break;
				}
				else if (indexDatasetArr[j][k] > -1)	sameF	= false;

			}
			if (k >= n && (!sameF || (sameF && i < j))) // if sameF is true, there are two occasions. we can take when i < j or i > j
				findExactClusterPart(i, j);
		}
	}

//	for (unsigned int i = 0, m = clusterExactPairArr.size(); i < m; ++i)
//		{
//			cout << i << ": " << (clusterExactPairArr[i].size()) << " : ";
//			for (unsigned int j = 0, n = clusterExactPairArr[i].size(); j < n; ++j)
//				cout << clusterExactPairArr[i][j] << "\t";
//			cout << endl;
//		}

	vector<int> tempArr;
	for (unsigned int i = 0, m = clusterExactPairArr.size(); i < m; ++i)
	{
		if (clusterExactPairArr[i][0] >= 0)
			for (unsigned int j = 1, n = clusterExactPairArr[i].size(); j < n; ++j)
			{
				if ((clusterExactPairArr[i][0] - clusterExactPairArr[i][j]) != 0)
				{
					appendExactInd(clusterExactPairArr[i][j], tempArr);
					for (unsigned int k = 0; k < tempArr.size(); ++k)
						if (tempArr[k] != clusterExactPairArr[i][0])
							clusterExactPairArr[i].push_back(tempArr[k]);
				}
			}
	}


	for (unsigned int i = 0, m = clusterExactPairArr.size(); i < m; ++i)
	{
		if (clusterExactPairArr[i][0] >= 0)
		{
			if (clusterExactPairArr[i].size() > 1)
			{
				vector<int> rowArr;
				for (unsigned int j = 0, n = clusterExactPairArr[i].size(); j < n; ++j)
				{
					if (find(rowArr.begin(), rowArr.end(), clusterExactPairArr[i][j]) - rowArr.begin() >= rowArr.size())
						rowArr.push_back(clusterExactPairArr[i][j]);
				}
				clusterExactIndArr.push_back(rowArr);
			}
			else
				clusterExactIndArr.push_back(clusterExactPairArr[i]);
		}
	}

//	for (unsigned int i = 0, m = clusterExactIndArr.size(); i < m; ++i)
//	{
//		cout << i << ": " << (clusterExactIndArr[i].size()) << " : ";
//		for (unsigned int j = 0, n = clusterExactIndArr[i].size(); j < n; ++j)
//			cout << clusterExactIndArr[i][j] << "\t";
//		cout << endl;
//	}
}

void appendExactInd(int pos, vector<int>& tempArr)
{
	if (clusterExactPairArr[pos][0] < 0)
		return;

	for (unsigned int i = 1; i < clusterExactPairArr[pos].size(); ++i)
		tempArr.push_back(clusterExactPairArr[pos][i]);

	clusterExactPairArr[pos][0]	= -1;
}

// find exact clusters for each pair
void findExactClusterPart(int indSub, int indSup)
{

	int indField;
	vector<StrPacket> tempArr;
	vector<int> truncArr;

	for (unsigned int i = 0, m = weightArr.size(); i < m; ++i)
		truncArr.push_back(-1);

	for (unsigned int i = 0, m = attrArr[TYPE_TRUNC].size(); i < m; ++i)
		truncArr[attrArr[TYPE_TRUNC][i][0]] = attrArr[TYPE_TRUNC][i][1];

	// strpacket creation for sub data set
	int indStartSub 	= recordStartIndArr[indSub];
	int indLastSub	 	= (indSub + 1) < fileNameArr.size()?  recordStartIndArr[indSub + 1] - 1 : recordArr.size() - 1;

	for (int i = indStartSub; i <= indLastSub; ++i)
	{
		string strCC	= ""; // concatenated string

		for (unsigned int j = 0, m = indexDatasetArr[indSub].size(); j < m; ++j)
		{
			indField	= indexDatasetArr[indSub][j];
			if (indField >= 0 && weightArr[j] > 0)
			{
				if (truncArr[j] > -1 && recordArr[i][indField].length() >= (unsigned)truncArr[j])
					strCC.append(recordArr[i][indField].substr(0, truncArr[j])).append("#");
				else
					strCC.append(recordArr[i][indField]).append("#");
			}
		}

		tempArr.push_back((StrPacket){strCC, i});
	}

	// strpacket creation for super data set
	int indStartSup 	= recordStartIndArr[indSup];
	int indLastSup	 	= (indSup + 1) < fileNameArr.size()?  recordStartIndArr[indSup + 1] - 1 : recordArr.size() - 1;

	int indexSupArr[indexDatasetArr[indSub].size()];
	for (int j = 0; j < indexDatasetArr[indSub].size(); ++j)
	{
		if (indexDatasetArr[indSub][j] > -1)
			indexSupArr[j]	= indexDatasetArr[indSup][j];
		else
			indexSupArr[j]	= -1;
	}

	for (int i = indStartSup; i <= indLastSup; ++i)
	{
		string strCC	= ""; // concatenated string

		for (unsigned int j = 0, m = indexDatasetArr[indSup].size(); j < m; ++j)
		{
			indField	= indexSupArr[j];
			if (indField >= 0 && weightArr[j] > 0)
			{
				if (truncArr[j] > -1 && recordArr[i][indField].length() >= (unsigned)truncArr[j])
					strCC.append(recordArr[i][indField].substr(0, truncArr[j])).append("#");
				else
					strCC.append(recordArr[i][indField]).append("#");
			}
		}

		tempArr.push_back((StrPacket){strCC, i});
	}

	cout << indSub << " in " << indSup << " record " << tempArr.size() << endl;
	//cout << indStartSub << ":" << indLastSub << "\t" << indStartSup << ":" << indLastSup << endl;

	vector<StrPacket> strSortedArr	= sortMSD(tempArr);
	//for (unsigned int k = 0; k < strSortedArr.size(); ++k)
		//	cout << strSortedArr[k].str << " index " << strSortedArr[k].ind << endl;

	vector<int> clusterRowArr;
	clusterRowArr.push_back(strSortedArr[0].ind);
	for (unsigned int i = 1, m = strSortedArr.size(); i < m; ++i)
	{
		if (strSortedArr[i].str.length() > 5 && strSortedArr[i].str == strSortedArr[i - 1].str)
		{
			clusterRowArr.push_back(strSortedArr[i].ind);
		}
		else
		{
			if (clusterRowArr.size() > 1)
			{
				int idSup	= -1;
				for (unsigned int j = 0; j < clusterRowArr.size(); ++j)
					if (clusterRowArr[j] < indStartSub || clusterRowArr[j] > indLastSub)
					{
						idSup	= clusterRowArr[j];
						break;
					}
				if (idSup > -1)
				{
					for (unsigned int j = 0; j < clusterRowArr.size(); ++j)
						if (clusterRowArr[j] >= indStartSub && clusterRowArr[j] <= indLastSub)
							clusterExactPairArr[idSup].push_back(clusterRowArr[j]);
				}
			}

			clusterRowArr.clear();
			clusterRowArr.push_back(strSortedArr[i].ind);
		}
	}

	if (clusterRowArr.size() > 1)
	{
		int idSup	= -1;
		for (unsigned int j = 0; j < clusterRowArr.size(); ++j)
			if (clusterRowArr[j] < indStartSub || clusterRowArr[j] > indLastSub)
			{
				idSup	= clusterRowArr[j];
				break;
			}
		if (idSup > -1)
		{
			for (unsigned int j = 0; j < clusterRowArr.size(); ++j)
			{
				if (clusterRowArr[j] >= indStartSub && clusterRowArr[j] <= indLastSub)
					clusterExactPairArr[idSup].push_back(clusterRowArr[j]);
			}
		}
	}
}



// find approximate cluster using approximate match by blocking, generating edgelist followed by findig connected components
// RURU
void findApproxCluster(vector<int>& rootArr)
{
	clock_t currTS2	= clock();
	
	
	vector<int> coverSetArr;

	for (int i = 0; i < fileNameArr.size(); ++i)
		coverSetArr.push_back(0);

	int count;

	for (int i = 0; i < blockFieldArr[0].size(); ++i)
	{
		count	= 0;
		for (int j = 0; j < indexDatasetArr.size(); ++j)
			if (indexDatasetArr[j][blockFieldArr[0][i]] >= 0)
			{
				++count;
				coverSetArr[j]	= 1;
			}

		if (count < 2)
			continue;
		
		vector<vector<int>> blockArr;
		createBlock(blockFieldArr[0][i], blockFieldArr[1][i] - NUMBER_LMER + 1, blockFieldArr[2][i], blockArr);
		createClusterEdgeList(blockArr);

		for (int j = 0; j < coverSetArr.size(); ++j)
			if (coverSetArr[j] < 1)
				break;

		//if (j >= coverSetArr.length)
						//break;

		//if (i == 0)
			//break;
	}
	
	cout << "TIME " << (double)(clock() - currTS2) / CLOCKS_PER_SEC << endl;
	findConnComp(rootArr); // find connected components on the graph generated by edgelist as edges and records as points
}



// create blocks of records using LMER(here LMER = 3) characters of last name
// Ruru

void createBlock(int indBlockField, int lmerUsed, int type, vector<vector<int>>& blockArr)
{
	cout << clusterExactIndArr.size() << " createBlock() : " << indBlockField << "\ttype " << type << " kmer: " << lmerUsed << endl;
	if (lmerUsed < 3)
		lmerUsed	= 3;
	int blockTotal 	= pow(ALPHABET_SIZE_LIST[type], lmerUsed);
	string strSample;

	if (type == 0) // alphabet is english alphabet
		strSample	= "aaaaaaaaaa"; // enough amount of characters for empty string (here 10)
	else
		strSample	= "0000000000";

	
	blockArr.resize(blockTotal);
	
	vector<string> record;
	int indFieldDataset, blockID;
	string blockFieldStr;

	int blkCount = 0;
	for (int i = 0; i < clusterExactIndArr.size(); ++i) {
		record	= recordArr[clusterExactIndArr[i][0]];
		// for(int k = 0; k< record.size(); k++) {
		// 	cout<< record[k]<<endl;
		// }
		indFieldDataset	= indexDatasetArr[atoi(record[record.size() - 1].c_str())][indBlockField];
		if (indFieldDataset < 0)
			continue;
		blockFieldStr	= record[indFieldDataset];
		int strLen	= blockFieldStr.length();
		//  cout << "block " << i << "\tX" << blockFieldStr << "X\t" << strLen << "\t" << record[1] << endl;
		if (strLen < lmerUsed) {
			blockFieldStr	= strSample.substr(0, lmerUsed - strLen) + blockFieldStr;
			strLen = lmerUsed;
		}
		
		//cout << i << "\t" << blockFieldStr << "\t" << strLen << endl;

		vector<int> codeRecordArr;
		//codeRecordArr.resize(strLen);

		for (int j = 0; j < blockFieldStr.length(); ++j)
		{
			if ( ( (type == 0) && !(((int) blockFieldStr[j] >= 97 && (int) blockFieldStr[j] <= 122)) )
			|| ( (type == 1) && !(((int) blockFieldStr[j] >= 48 && (int) blockFieldStr[j] <= 57)) )
			|| ( (type == 2) && (!((((int) blockFieldStr[j] >= 97 && (int) blockFieldStr[j] <= 122)) || (((int) blockFieldStr[j] >= 48 && (int) blockFieldStr[j] <= 57)))) )
			)
			{
				blockFieldStr	= blockFieldStr.substr(0, j) + blockFieldStr.substr(j + 1);
				--j;
			}
			else
				codeRecordArr.push_back((int) blockFieldStr[j]);
		}

		blockID=0;
		for (int j = 0; j < blockFieldStr.length() - lmerUsed + 1; ++j)
		{
			blockID	= 0;
			for (int k = 0; k < lmerUsed; ++k)
				if (type == 0)
						blockID	+= (codeRecordArr[j + k] - 97) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1);
				else if(type == 1)
					blockID	+= (codeRecordArr[j + k] - 48) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1);
				else if(type == 2)
				{
					if(codeRecordArr[j + k] >= 97)
						blockID	+= (codeRecordArr[j + k] - 97) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1);
					else
						blockID	+= (codeRecordArr[j + k] - 22) * (int) pow(ALPHABET_SIZE_LIST[type], lmerUsed - k - 1); // 48 - 26
				}

			//if (atoi(record[1].c_str()) == 242299882 || atoi(record[1].c_str()) == 242783653)
				//cout << i << ":" << blockFieldStr << " id:" << blockID << " len:" << blockFieldStr.length() << " " << j << "\t" << (blockFieldStr.length() - lmerUsed + 1) << "::" << record[record.size() - 1] << endl;
			if(!(blockID < 0 || blockID >= blockTotal)) {
				if (is_rand==true) {
					// RURU found it
					int rand_int;
					rand_int = rand() % 100;
					if (rand_int < rand_theshold) {
						blockArr[blockID].push_back(i);
					}
				} else {
					blockArr[blockID].push_back(i);
				}
				
				//cout<<"BlockID: "<< blockID << " cluster i: "<< i <<endl;
			}
				
		}

		blkCount	+= (blockFieldStr.length() - lmerUsed + 1);
	}
	cout << blockTotal << "blk count: " << blkCount << endl;
}

int c1 = 0, c2 = 0;
void createClusterEdgeList(vector<vector<int>>& blockArr)
{
	cout << "createClusterEdgeList" << endl;
	int	blockTotal	= blockArr.size();
	testVal	= 0;
	for (int i = 0; i < blockTotal; ++i)
	{
		if (blockArr[i].size() > 0)
		{
			//cout << i << " size " << blockArr[i].size() << endl;
			generateEdgilist(blockArr[i]);
		}
	}

	cout <<  c1 << " " << c2 << " link total " << testVal << " blk: " << blockTotal << endl;
}


 // generate edge list within a block

void generateEdgilist(vector<int>& blockRowArr)
{
	int blockItemTotal	= blockRowArr.size();

	vector<vector<string> > dataArr(blockItemTotal); // to make cache-efficient, keep records in a row
	for (int i = 0; i < blockItemTotal; ++i)
		dataArr[i]	= recordArr[clusterExactIndArr[blockRowArr[i]][0]];
	vector<int> vectorArr;
	vectorArr.resize(blockItemTotal, 0);

	int temp, posTemp;
	int n = 0;
	for (int i = 0; i < blockItemTotal; i++)
	{
		if (vectorArr[i] == 0)
		{c2++;
			vector<int> tempArr;
			tempArr.push_back(i);
			posTemp			= 0;
			vectorArr[i]	= (++n);
			while (posTemp < tempArr.size())
			{
				temp		= tempArr[posTemp++];
				for (int j = 0; j < blockItemTotal; j++)
				{
					if (vectorArr[j] == 0 && j != temp)
					{	
					c1++;
						if (isLinkageOk(dataArr[temp], dataArr[j]))
						{
							tempArr.push_back(j);
							vectorArr[j]	= n;
							vector<int> edge(2);
							edge[0]	= blockRowArr[temp], edge[1] = blockRowArr[j];
							edgeArr.push_back(edge);
						}
					}
				}
			}
		}
	}

	dataArr.clear();
}


//  find clusters as connected components in a graph where edges are connection among records and vertices are records

void findConnComp(vector<int>& parentArr)
{
	//cout << "findConnComp() " << edgeArr.size() << endl;

	//startCompT	= clock();
	int i, rootU, rootV, edgeTotal, pointTotal;
	vector<int> weightArr;

	pointTotal	= clusterExactIndArr.size();
	for(i = 0; i < pointTotal; ++i)
	{
		parentArr.push_back(i);
		weightArr.push_back(0);
	}

	edgeTotal	= edgeArr.size();
	//cout << edgeTotal <<  " et" << endl;
	for(i = 0; i < edgeTotal; ++i)
	{
		//cout << edgeArr[i][0] << ":" << endl;
		rootU	= findRoot(edgeArr[i][0], parentArr);
		rootV	= findRoot(edgeArr[i][1], parentArr);


		if(rootU != rootV)
		{
			//cout << i << " e " << rootU << ":" << rootV << " edge " << edgeArr[i][0] << ":" << edgeArr[i][1] << endl;
			makeEquivalent(rootU, rootV, parentArr, weightArr);
		}
	}
	cout << "findConnComp()" << endl;
}


//  combine exact clusters and approximate clusters

void findFinalCluster(vector<int>& rootArr)
{
	cout << "findFinalCluster() " << rootArr.size() << endl;

	int indCluster;

	clusterIndArr.assign(rootArr.size(), vector<int>());

	for (unsigned int i = 0; i < rootArr.size(); ++i)
	{
		indCluster	= findRoot(i, rootArr);
		//cout << i << " c " << indCluster << " s " << clusterIndArr.size() << endl;
		for (unsigned int j = 0; j < clusterExactIndArr[i].size(); ++j)
			clusterIndArr[indCluster].push_back(clusterExactIndArr[i][j]);
	}

	clusterExactIndArr.clear();
	rootArr.clear();

	int k = 0;
	for (unsigned int i = 0; i < clusterIndArr.size(); ++i)
	{
		if (clusterIndArr[i].size() > 0)
			k++;
	}

	cout << k << " cluster total" << endl;


//	int k = 0;
//	for (unsigned int i = 0; i < 2; ++i)
//	{
//		if (clusterIndArr[i].size() > 0)
//			cout << endl << "cluster " << (k++) << ":";
//		for (unsigned int j = 0; j < clusterIndArr[i].size(); ++j)
//			cout << "\t" << clusterIndArr[i][j];
//	}

}


void clusterGrp(vector<int>& recordIndSingleGrpArr, vector<Cluster>& clusterSingleArr)
{
	vector<vector<string> > recordGrpArr(recordIndSingleGrpArr.size());

	for (unsigned int i = 0; i < recordIndSingleGrpArr.size(); ++i)
	{
		recordGrpArr[i]	= recordArr[recordIndSingleGrpArr[i]];
		//if(recordArr.get(recordIndSingleGrpArr.get(i)).get(2).equalsIgnoreCase("gaines"))
			//System.out.println(recordIndSingleGrpArr.size() + " ind " + i + " ERSR " + recordArr.get(recordIndSingleGrpArr.get(i)));
	}

	int recordGrpTotal	= recordGrpArr.size();
	//cout << " grp " << recordGrpTotal << endl;;

	for (int i = 0; i < recordGrpTotal; i++)
	{
		Cluster cluster;
		cluster.initCluster(0, recordGrpArr[i]);
		clusterSingleArr.push_back(cluster);
	}



	vector<vector<int> > matArr, vecArr;
	generateMatrix(recordGrpArr, matArr);
/*
	for (unsigned int i = 0; i < matArr.size(); i++)
	{
		cout << "\n" << i << ":";
		for (unsigned int j = 0; j < matArr[i].size(); j++)
		{
			cout << matArr[i][j] << "\t";
		}
	}
	vecArr	= generateVector(matArr);
	for (unsigned int i = 0; i < vecArr.size(); i++)
	{
		cout << "\n" << i << ":";
		for (unsigned int j = 0; j < vecArr[i].size(); j++)
		{
			cout << vecArr[i][j] << "\t";
		}
	}
*/
	generateVector(matArr, vecArr);

	//cout << "Complete GenVec\n";
	int distCluster, distTemp;
	int idCluster1, idCluster2, idMinCluster, idMaxCluster;

	while (matArr.size() >= 2)
	{
		distCluster	= vecArr[0][0];
		idCluster1	= (int) vecArr[0][1];
		idCluster2	= (int) vecArr[0][2];
		//cout << distCluster << ":" << idCluster1 << "\t" << idCluster2 << endl;
		for (unsigned int i = 1; i < vecArr.size(); i++)
		{
			distTemp	= vecArr[i][0];
			if (distCluster > distTemp)
			{
				distCluster	= distTemp;
				idCluster1	= (int) vecArr[i][1];
				idCluster2	= (int) vecArr[i][2];
			}
		}
		//cout << distCluster << " ch:" << idCluster1 << "\t" << idCluster2 << endl;
		if (distCluster <= threshold)
		{
			if (idCluster1 > idCluster2)
			{
				idMinCluster	= idCluster2;
				idMaxCluster	= idCluster1;
			}
			else
			{
				idMinCluster	= idCluster1;
				idMaxCluster	= idCluster2;
			}
			//cout << idMaxCluster << "X" << idMinCluster << endl;
			Cluster newC;
			newC.initCluster(distCluster, clusterSingleArr[idMinCluster], clusterSingleArr[idMaxCluster]);
			clusterSingleArr.erase(clusterSingleArr.begin() + idMaxCluster);
			clusterSingleArr.erase(clusterSingleArr.begin() + idMinCluster);
			clusterSingleArr.push_back(newC);

			updateMatVec(matArr, vecArr, idMinCluster, idMaxCluster);
			//cout << matArr.size() << " size " << vecArr.size() << " id " << idMinCluster << ":" << idMaxCluster << endl;

			/*
			for (unsigned int i = 0; i < vecArr.size(); i++)
				{
					cout << "\n" << i << ":";
					for (unsigned int j = 0; j < vecArr[i].size(); j++)
					{
						cout << vecArr[i][j] << "\t";
					}
				}
			cout << endl;
			*/
			//break;
		}
		else
			break;
	}

	//cout << "finding clusters completed " << clusterSingleArr.size() << endl;

	vector<string> recordItem;
	int fileInd, k;


	//System.out.println("finding clusters completed " + indGrp + " size: " + clusterSingleArr.size());

	//clusterSingleArr.get(0).PrintCluster();
	if (clusterSingleArr.size() > 1)
		postprocessCluster(clusterSingleArr);
	//cout << "post finding clusters completed " << clusterSingleArr.size() << endl;
	//return clusterSingleArr;

	for (unsigned int i = 0; i < clusterSingleArr.size(); ++i)
	{
		if (clusterSingleArr[i].itemArr.size() > 1)
			for (unsigned int j = 0; j < clusterSingleArr[i].itemArr.size(); ++j)
			{
				recordItem	= clusterSingleArr[i].itemArr[j];
				fileInd		= atoi(recordItem[recordItem.size() - 1].c_str());

				for (k = 0; k < priorityFieldArr.size(); ++k)
				{
//					if(fileInd == 3 && recordItem.at(2) == "mccall")
//						cout << priorityFieldArr[k] << " p " << k << " ind " << indexDatasetArr[fileInd][priorityFieldArr[k]] << "\t" << priorityFieldArr.size() <<   endl;

					if ((indexDatasetArr[fileInd][priorityFieldArr[k]] >= 0) && (recordItem[indexDatasetArr[fileInd][priorityFieldArr[k]]].length() > 0))
						break;
				}

				if (k >= priorityFieldArr.size())
				{
					recordItem	= clusterSingleArr[i].itemArr[j];
					//cout << recordItem[0] << " record " << recordItem[2] << " file: " << fileInd << endl;
					clusterSingleArr[i].itemArr.erase(clusterSingleArr[i].itemArr.begin() + j);
					//if(recordItem.get(1).equalsIgnoreCase("shaeleigh"))
						//System.out.println("ERROR " + recordItem);
					Cluster newC;
					newC.initCluster(0, recordItem);
					clusterSingleArr.push_back(newC);
				}
			}
	}
}



// generate 2d matrix of distances among records in recordGrpArr
void generateMatrix(vector<vector<string>>& recordGrpArr, vector<vector<int> >& matArr)
{
	;
	int recordGrpTotal;

	recordGrpTotal	= recordGrpArr.size();

	for (int i = 0; i < recordGrpTotal; i++)
	{
		vector<int> rowMatArr;

		for (int j = 0; j < i; j++)
		{
			rowMatArr.push_back(linkage(recordGrpArr[i], recordGrpArr[j]));
		}
		matArr.push_back(rowMatArr);
	}
}


// generate vector of entries having min distance for each index
void generateVector(vector<vector<int> >& matArr, vector<vector<int> >& vecArr)
{
	int distMin, indThis, indOther;

	for (unsigned int i = 0; i < matArr.size(); i++)
	{
		indThis 	= i;
		if (i != (matArr.size()-1))
		{
			distMin		= matArr[i+1][i];

			indOther	= i + 1;
		}
		else
		{
			distMin 	= matArr[i][0];
			//indThis		= i;
			indOther	= 0;
		}

		for (unsigned int j = 0; j < matArr.size(); j++)
		{
			if (i > j)
			{
				if (distMin > matArr[i][j])
				{
					distMin		= matArr[i][j];
					//indThis		= i;
					indOther	= j;
				}
			}
			else if (i < j)
			{
				if (distMin > matArr[j][i])
				{
					distMin		= matArr[j][i];
					//indThis		= i;
					indOther	= j;
				}
			}
		}

		vector<int> obj;
		obj.push_back(distMin);
		obj.push_back(indThis);
		obj.push_back(indOther);
		vecArr.push_back(obj);
	}
}


// update matArr and vecArr after merging two clusters
void updateMatVec(vector<vector<int> > &matArr, vector<vector<int> > &vecArr, int indCluster1, int indCluster2)
{
	//cout << "updateMatVec\n";
	int indMin, indMax;
	int dist1, dist2;
	vector<int> rowMatArr;

	indMin	= indCluster1;
	indMax	= indCluster2;

	for (unsigned int i = 0; i < matArr.size(); i++)
	{
		if(i != indCluster1 && i != indCluster2)
		{
			if (i > indCluster1)
				dist1	= matArr[i][indCluster1];
			else
				dist1	= matArr[indCluster1][i];

			if (i > indCluster2)
				dist2	= matArr[i][indCluster2];
			else
				dist2	= matArr[indCluster2][i];

			rowMatArr.push_back(max(dist1, dist2));
		}
	}

	//remove rows and cols from m
	matArr.erase(matArr.begin() + indMax);
	matArr.erase(matArr.begin() + indMin);

	for (unsigned int i = 0; i < matArr.size(); i++)
	{
		if (matArr[i].size() >= indMax)
			matArr[i].erase(matArr[i].begin() + indMax);
		if (matArr[i].size() >= indMin)
			matArr[i].erase(matArr[i].begin() + indMin);
	}

	//add new row for r1+r2 to m
	matArr.push_back(rowMatArr);

	//get min distance and min distance index for r1+r2

	if (rowMatArr.size() != 0)
	{
		int distMinNew;
		int idThisNew, idOtherNew, idTemp;


		distMinNew	= rowMatArr[0];
		idOtherNew	= 0;
		for (unsigned int i = 0; i < rowMatArr.size(); i++)
		{
			if (distMinNew > rowMatArr[i])
			{
				distMinNew	= rowMatArr[i];
				idOtherNew	= i;
			}
		}

		//remove r1 and r2 from v
		vecArr.erase(vecArr.begin() + indMax);
		vecArr.erase(vecArr.begin() + indMin);

		//add new info for r1+r2
		idThisNew	= vecArr.size(); //insert location for info of r1+r2
		//static const float obj[] = {distMinNew, (float)idThisNew, (float)idOtherNew};
		vector<int> rowVec;
		rowVec.push_back(distMinNew);
		rowVec.push_back(idThisNew);
		rowVec.push_back(idOtherNew);

		//cout << indMax << " mm " << indMin << "\t" << distMinNew << ":" << idThisNew << " " << idOtherNew << endl;
		 //cout << rowVec[0] << ":" << rowVec[1] << " " << rowVec[2] << endl;

		for (unsigned int i = 0; i < vecArr.size(); i++)
		{
			idTemp	= vecArr[i][1];
			if (idTemp > indMax)
				vecArr[i][1]	= (idTemp - 2); //update v[i][1] since remove of r1 and r2
			else if (idTemp > indMin)
				vecArr[i][1]	= (idTemp - 1);

			idTemp	= (int) vecArr[i][2];
			//System.out.println(i + "->" + v.get(i).get(1) + ":" + v.get(i).get(2));
			if (idTemp == indCluster1 || idTemp == indCluster2)
			{
				vecArr[i][0]	= matArr[matArr.size() - 1][round(vecArr[i][1])];
				vecArr[i][2]	= idThisNew; // update nearest neighbor
			}
			else
			{
				if (idTemp > indMax)
					vecArr[i][2]	= (idTemp - 2); //update v[i][2] since remove of r1 and r2
				else if (idTemp > indMin)
					vecArr[i][2]	= (idTemp - 1);

			}
			//System.out.println(i + "->" + v.get(i).get(1) + ":" + v.get(i).get(2));
		}


		vecArr.push_back(rowVec);
	}
}



void postprocessCluster(vector<Cluster>& clusterSingleArr)
{

	//System.out.println("postprocessCluster: " + clusterSingleArr.size());

	//cout << "postprocessCluster: " << clusterSingleArr.size() << endl;

	int l, maxVal, maxInd;
	Cluster cluster1, cluster2;
	vector<string> record;
	vector<vector<int> > usedAttrArr(clusterSingleArr.size());
	vector<int> isInClusterArr;

	for (unsigned int i = 0; i < clusterSingleArr.size(); ++i)
	{
		for (unsigned int j = 0; j < weightArr.size(); ++j)
			usedAttrArr[i].push_back(0);

		cluster1	= clusterSingleArr[i];
		for (unsigned int j = 0; j < cluster1.itemArr.size(); ++j)
		{
			record	= cluster1.itemArr[j];
			l		= atoi(record[record.size() - 1].c_str());

			for (unsigned int k = 0; k < indexDatasetArr[l].size(); ++k)
			{
				if (indexDatasetArr[l][k] > -1)
					usedAttrArr[i][k]	= 1;
			}
		}
	}

	//cout << "pos \n";
	for (unsigned int i = 0; i < clusterSingleArr.size(); ++i)
	{

		//cluster1	= clusterSingleArr[i];

		for (unsigned int j = 0; j < clusterSingleArr[i].itemArr.size(); )
		{
			//cout << "bpos " << i << ":" << j << endl;
			isInClusterArr.assign(clusterSingleArr.size(), 0);
			isInClusterArr[i]	= 1;

			record	= clusterSingleArr[i].itemArr[j];

//			if(record[1]== "242605281" || record[1] == "242603923")
//			{
//				cout << i << " " << record[0] << ":" << record[1] << " size " << clusterSingleArr.size() << "::" << clusterSingleArr[i].itemArr.size() << endl;
//				for (unsigned int v = 0; v < clusterSingleArr.size(); v++)
//				{
//					cout << (clusterSingleArr[v].itemArr[0])[0] << endl;
//				}
//			}

			for (unsigned int k = 0; k < clusterSingleArr.size(); ++k)
			{
				if(i == k)
					continue;

				cluster2	= clusterSingleArr[k];

				if(isInCluster(cluster2, record)) // check if this record also is in cluster2
				{
					//if(record[1]== "242605281" || record[1] == "242603923")
						//cout << i << " in " << k << " -- " << record[2];
					isInClusterArr[k]	= 1;
				}
			}

			//cout << "pos " << i << ":" << j << endl;
			vector<vector<int> > usedAttrThisArr;
			processUsedAttr(usedAttrArr, isInClusterArr, atoi(record[record.size() - 1].c_str()), usedAttrThisArr);
			//cout << "apos " << i << endl;
			int fileIndThis	= atoi(record[record.size() - 1].c_str());
			for (unsigned int l = 0; l < clusterSingleArr.size(); ++l)
			{
				if(isInClusterArr[l] > 0)
				{
					for (int t = 0; t < clusterSingleArr[l].itemArr.size(); ++t)
					{
						vector<string> recordTemp	= clusterSingleArr[l].itemArr[t];
						int fileIndTemp	= atoi(recordTemp[recordTemp.size() - 1].c_str());
						for (unsigned int k = 0; k < priorityFieldArr.size(); ++k)
						{
							if (indexDatasetArr[fileIndThis][priorityFieldArr[k]] >= 0 && indexDatasetArr[fileIndTemp][priorityFieldArr[k]] >= 0)
							{
								if (record[indexDatasetArr[fileIndThis][priorityFieldArr[k]]] == recordTemp[indexDatasetArr[fileIndTemp][priorityFieldArr[k]]])
									isInClusterArr[l]	+= (priorityFieldArr.size() - k) * 100;
								else
									isInClusterArr[l]	-= (priorityFieldArr.size() - k) * 100;
							}
						}
					}
				}
			}

			maxVal	= isInClusterArr[i];
			maxInd	= i;
			for (unsigned int k = 0; k < isInClusterArr.size(); ++k)
			{
				if ((i != k) && (isInClusterArr[k] >= maxVal))
				{
					maxVal	= isInClusterArr[k];
					maxInd	= k;
				}
			}

			//cout << "Xpos " << i << ":" << j << " max " << maxInd << " size " << clusterSingleArr[i].itemArr.size() <<  endl;

			if(maxInd != i)
			{
				//if(record.get(2).equals("gaines"))
					//System.out.println(i + " to " + maxInd + " rr " + record);
				record	= clusterSingleArr[i].itemArr[j];
				clusterSingleArr[i].itemArr.erase(clusterSingleArr[i].itemArr.begin() + j);
				clusterSingleArr[maxInd].itemArr.push_back(record);

			}
			else
				j++;
		}
	}
}

bool isInCluster(Cluster& cluster, vector<string>& record)
{
	//System.out.println("isInCluster: ");

	unsigned int i;
	for (i = 0; i < cluster.itemArr.size(); ++i)
		if (!isLinkageOk(record, cluster.itemArr[i]))
			break;

	if (i < cluster.itemArr.size())
		return false;
	else
		return true;
}

void processUsedAttr(vector<vector<int> >& usedAttrArr, vector<int>& isInClusterArr, int indDataset, vector<vector<int> >& usedThisAttrArr)
{
	//cout << "processUsedAttr: " << endl;

	usedThisAttrArr.resize(usedAttrArr.size());

	for (unsigned int i = 0; i < usedAttrArr.size(); ++i)
	{
		usedThisAttrArr[i].assign(weightArr.size(), 0);

		if (isInClusterArr[i] != 0)
		{
			for (unsigned int j = 0; j < weightArr.size(); ++j)
			{
				if (indexDatasetArr[indDataset][j] >= 0)
					usedThisAttrArr[i][j]	= (usedAttrArr[i][j]);
			}
		}
	}
}



void output(int outputType)
{
	cout << "\noutput()" << endl;

	string fileNameStr	= outDir;
	if (outputType == OUTPUT_FINAL)
		fileNameStr.append("OutFinal.txt");
	else
		fileNameStr.append("OutSingle.txt");
	cout << "output " << fileNameStr << endl;



	unsigned int clusterInd;
	vector<string> recordRow;
	ofstream outFile;
	outFile.open(fileNameStr.c_str(), ofstream::out);

	vector<vector<int> > typeClusterIndArr(4);
	vector<int> typeClusterArr;
	typeClusterArr.resize(4, 0);
	clusterInd	= 0;
	if (outputType == OUTPUT_FINAL)
	{

		cout << "Final Output: " << endl;
		for (unsigned int i = 0; i < clusterArr.size(); ++i)
		{
			if (clusterArr[i].itemArr.size() <= 0)
				continue;
			vector<int> countArr, indArr;

			++clusterInd;
			outFile << "cluster " << clusterInd << ": " << endl;
			for (unsigned int j = 0; j < clusterArr[i].itemArr.size(); ++j)
			{
				recordRow	= clusterArr[i].itemArr[j];
				int recordID	= atoi(recordRow[0].c_str());
				int pos		= find(indArr.begin(), indArr.end(), recordID) - indArr.begin();
				if (pos < countArr.size())
					countArr[pos]++;
				else
				{
					indArr.push_back(recordID);
					countArr.push_back(1);
				}
				for (unsigned int k = 0; k < recordRow.size(); ++k)
				{
					outFile << recordRow[k] << "\t";
				}
				outFile << endl;


			}
			outFile << endl;

			if (indArr.size() == 1) // T1 ot T2
			{
				if (countArr[0] == fileNameArr.size())
					typeClusterArr[0]	+= clusterArr[i].itemArr.size();
				else
				{
					typeClusterIndArr[1].push_back(clusterInd);
					typeClusterArr[1]	+= clusterArr[i].itemArr.size();
				}
			}
			else
			{
				int t;
				for (t = 0; t < countArr.size(); ++t)
					if (countArr[t] == fileNameArr.size())
						break;
				if (t < countArr.size())
				{
					typeClusterIndArr[2].push_back(clusterInd);
					typeClusterArr[2]	+= clusterArr[i].itemArr.size();
				}
				else
				{
					typeClusterIndArr[3].push_back(clusterInd);
					typeClusterArr[3]	+= clusterArr[i].itemArr.size();
				}
			}
		}
	}
	else
	{
		cout << "SingleLinkageOutput: " << endl;
		for (unsigned int i = 0; i < clusterIndArr.size(); ++i)
		{
			if (clusterIndArr[i].size() <= 0)
				continue;
			vector<int> countArr, indArr;

			++clusterInd;
			outFile << "cluster " << clusterInd << ": " << endl;
			for (unsigned int j = 0; j < clusterIndArr[i].size(); ++j)
			{
				recordRow	= recordArr[clusterIndArr[i][j]];
				int recordID	= atoi(recordRow[0].c_str());
				int pos		= find(indArr.begin(), indArr.end(), recordID) - indArr.begin();
				if (pos < countArr.size())
					countArr[pos]++;
				else
				{
					indArr.push_back(recordID);
					countArr.push_back(1);
				}
				for (unsigned int k = 0; k < recordRow.size(); ++k)
				{
					outFile << recordRow[k] << "\t";
				}
				outFile << endl;


			}
			outFile << endl;

			if (indArr.size() == 1) // T1 ot T2
			{
				if (countArr[0] == fileNameArr.size())
					typeClusterArr[0]	+= clusterIndArr[i].size();
				else
				{
					typeClusterIndArr[1].push_back(clusterInd);
					typeClusterArr[1]	+= clusterIndArr[i].size();
				}
			}
			else
			{
				int t;
				for (t = 0; t < countArr.size(); ++t)
					if (countArr[t] == fileNameArr.size())
						break;
				if (t < countArr.size())
				{
					typeClusterIndArr[2].push_back(clusterInd);
					typeClusterArr[2]	+= clusterIndArr[i].size();
				}
				else
				{
					typeClusterIndArr[3].push_back(clusterInd);
					typeClusterArr[3]	+= clusterIndArr[i].size();
				}
			}
		}
	}

	outFile  << endl << "Total Cluster: " << clusterInd << endl;
	outFile.close();

	int typeCountSum	= 0;
	int recordPerFile	= floor(recordTotalCMD / fileNameArr.size());
	for (int t = 0; t < 4; ++t)
	{
		typeCountSum	+= typeClusterArr[t];
		if (typeClusterArr[t] > 0)
		{
			cout << "type " << (t + 1) << ": " << typeClusterArr[t] << " (%total: " << ((float)typeClusterArr[t] / recordTotal * 100) << ")" << endl;
//			if (t > 0)
//			{
//			for (int i = 0; i < typeClusterIndArr[t].size(); ++i)
//				cout << typeClusterIndArr[t][i] << "\t";
//			cout << endl;
//			}
		}
		else
			cout << "type " << (t + 1) << ": " << typeClusterArr[t] << " (%total: " << 0 << ")" << endl;

	}

	cout  << "Total Cluster: " << clusterInd << "\tType Count: " << typeCountSum << "\n\n";
}


bool isLinkageOk(vector<string>& a, vector<string>& b)
{
	int w	= calculateDistAll(a, b);
	return (threshold - w >= 0)? true : false;
}

int linkage(vector<string>& a, vector<string>& b)
{
	int w	= calculateDistAll(a, b);
	return w;
}

int calculateDistAll(vector<string>& a, vector<string>& b)
{
	int w	= 0;
	//cout << "linkage " << len << endl;
	for (int i = 0; i < 1; ++i)
	{
		if(i == TYPE_EDIT)
			w	+= calculateEditDist(a, b, attrArr[TYPE_EDIT], threshold - w);
		else if(i == TYPE_REVERSAL)
			w	+= calculateRevDist(a, b, attrArr[TYPE_REVERSAL], threshold - w);
		else if(i == TYPE_TRUNC)
			w	+= calculateTruncDist(a, b, attrArr[TYPE_TRUNC], threshold - w);
	}
	return w;
}

int compareAttLen(vector<string>& a, vector<string>& b)//length of compared attributes of a and b records
{
	int i_len = 0;
	int j_len = 0;

	for (int k = 0; k < attrArr[TYPE_EDIT].size(); k++)
	{
		int i_f = atoi(a[a.size()-1].c_str());
		int j_f = atoi(b[b.size()-1].c_str());

		int ind = attrArr[TYPE_EDIT][k][0];

		int i_index = indexDatasetArr[i_f][ind];
		int j_index = indexDatasetArr[j_f][ind];

		if (i_index == -1 || j_index == -1 || a[i_index].length() == 0 || b[j_index].length() == 0) continue;

		//System.out.println(a + ":" + i_index + " file: " + a.size());
		i_len += a[i_index].length();
		j_len += b[j_index].length();

	}

	for (int k = 0; k < attrArr[TYPE_REVERSAL].size(); k++)
	{
		int i_f = atoi(a[a.size()-1].c_str());
		int j_f = atoi(b[b.size()-1].c_str());

		int ind = attrArr[TYPE_REVERSAL][k][0];

		int i_index = indexDatasetArr[i_f][ind];
		int j_index = indexDatasetArr[j_f][ind];

		if (i_index == -1 || j_index == -1) continue;

		i_len += a[i_index].length();
		j_len += b[j_index].length();


		ind = attrArr[TYPE_REVERSAL][k][1];

		i_index = indexDatasetArr[i_f][ind];
		j_index = indexDatasetArr[j_f][ind];

		if (i_index == -1 || j_index == -1 || a[i_index].length() == 0 || b[j_index].length() == 0) continue;

		i_len += a[i_index].length();
		j_len += b[j_index].length();
	}

	for (int k = 0; k < attrArr[TYPE_TRUNC].size(); k++)
	{
		int i_f = atoi(a[a.size()-1].c_str());
		int j_f = atoi(b[b.size()-1].c_str());

		int ind = attrArr[TYPE_TRUNC][k][0];

		int i_index = indexDatasetArr[i_f][ind];
		int j_index = indexDatasetArr[j_f][ind];

		if (i_index == -1 || j_index == -1 || a[i_index].length() == 0 || b[j_index].length() == 0) continue;

		i_len += a[i_index].length();
		j_len += b[j_index].length();

	}

	return max(i_len, j_len);
}


int calculateEditDist(vector<string>& a, vector<string>& b, vector<vector<int> >& compareAtt, int threshRem)
{
	//cout << "calculateEditDist " << compareAtt.size() << endl;
	int w, ind, temp, a_ind, b_ind;
	string s1, s2;

	w 	= 0;

	int al_ind = atoi(a[a.size() - 1].c_str()), bl_ind = atoi(b[b.size() - 1].c_str());
	for (int g = 0; g < compareAtt.size(); g++)
	{
		ind 	= compareAtt[g][0];
		a_ind	= indexDatasetArr[al_ind][ind];
		b_ind	= indexDatasetArr[bl_ind][ind];

		if (a_ind == -1 || b_ind == -1 || a[a_ind].length() == 0 || b[b_ind].length() == 0)
			continue;

		c2++;
		s1 		= a[a_ind];
		s2 		= b[b_ind];
		temp 	= calculateBasicED(s1, s2, threshRem) * weightArr[ind];
		//if(a.get(2).equals("mercer") && b.get(2).equals("mercer"))
			//System.out.println(s1 + ":::" + s2 + " t " + temp);
		w 			+= temp;
		threshRem 	= threshRem-temp;
	}

	return w;
}


int calculateRevDist(vector<string>& a, vector<string>& b, vector<vector<int> >& compareAtt, int threshRem)
{
	//cout << "calculateRevDist " << compareAtt.size() << endl;

	int i, w, w1, w2, ind, ind1, ind2;
	string str1, str2;

	w	= 0;
	w1	= 0;
	w2	= 0;

	for(i = 0; i < compareAtt.size(); ++i)
	{
		ind		= compareAtt[i][0];
		ind1	= indexDatasetArr[atoi(a[a.size() - 1].c_str())][ind];
		ind2	= indexDatasetArr[atoi(b[b.size() - 1].c_str())][ind];
		if (ind1 == -1 || ind2 == -1 || a[ind1].length() == 0 || b[ind2].length() == 0)
			continue;
		str1	= a[ind1];
		str2	= b[ind2];
		w1		= calculateBasicED(str1, str2, threshRem - w);

		ind		= compareAtt[i][1];
		ind1	= indexDatasetArr[atoi(a[a.size() - 1].c_str())][ind];
		ind2	= indexDatasetArr[atoi(b[b.size() - 1].c_str())][ind];
		if (ind1 == -1 || ind2 == -1 || a[ind1].length() == 0 || b[ind2].length() == 0)
			continue;
		str1	= a[ind1];
		str2	= b[ind2];
		w1		= calculateBasicED(str1, str2, threshRem - w);

		ind		= compareAtt[i][0];
		ind1	= indexDatasetArr[atoi(a[a.size() - 1].c_str())][ind];
		ind		= compareAtt[i][1];
		ind2	= indexDatasetArr[atoi(b[b.size() - 1].c_str())][ind];
		if (ind1 == -1 || ind2 == -1 || a[ind1].length() == 0 || b[ind2].length() == 0)
			continue;
		str1	= a[ind1];
		str2	= b[ind2];
		w2		= calculateBasicED(str1, str2, threshRem - w);

		ind		= compareAtt[i][1];
		ind1	= indexDatasetArr[atoi(a[a.size() - 1].c_str())][ind];
		ind		= compareAtt[i][0];
		ind2	= indexDatasetArr[atoi(b[b.size() - 1].c_str())][ind];
		if (ind1 == -1 || ind2 == -1 || a[ind1].length() == 0 || b[ind2].length() == 0)
			continue;
		str1	= a[ind1];
		str2	= b[ind2];
		w2		+= calculateBasicED(str1, str2, threshRem - w);

		w		+= min(w1, w2);
	}

	return w;
}

int calculateTruncDist(vector<string>& a, vector<string>& b, vector<vector<int> >& compareAtt, int threshRem)
{
	//cout << "calculateTruncDist " << compareAtt.size() << endl;

	int w, indA, indB, indSetA, indSetB, count, ind;
	string str1, str2;

	w		= 0;
	indSetA	= atoi(a[a.size() - 1].c_str());
	indSetB	= atoi(b[b.size() - 1].c_str());

	for (unsigned int i = 0; i < compareAtt.size(); ++i)
	{
		ind		= compareAtt[i][0];
		count	= compareAtt[i][1];

		indA	= indexDatasetArr[indSetA][ind];
		if (indA == -1 || a[indA].length() == 0)
			continue;

		indB	= indexDatasetArr[indSetB][ind];
		if (indB == -1 || b[indB].length() == 0)
			continue;

		str1	= a[indA];
		str2	= b[indB];
		if (str1.length() >= count && str2.length() >= count)
		{
			str1	= str1.substr(0, count);
			str2	= str2.substr(0, count);
		}
		else
		{
			int len		= min(str1.length(), str2.length());
			str1	= str1.substr(0, len);
			str2	= str2.substr(0, len);
		}

		w	+= calculateBasicED(str1, str2, threshRem - w);
	}

	return w;
}



int calculateBasicED(string& str1, string& str2, int threshRem)
{
	int dist	= threshRem;

	if(abs((int)(str1.length() - str2.length())) > dist)
		return threshold + 1;
	else if(str1.compare(str2) == 0)
		return 0;
	else if(dist == 0)
		return threshold + 1;
	else if((2 * dist + 1) >= max(str1.length(), str2.length()))
		return calculateBasicED2(str1, str2, dist);
	else
	{
		string s1, s2;
		int row, col, diagonal;
		int i, j;
		vector<vector<int> > matArr;

		if (str1.length() > str2.length())
		{
			s1 = str2;
			s2 = str1;
		}
		else
		{
			s1 = str1;
			s2 = str2;
		}

		row	 		= s1.length() + 1;
		col 		= 2 * dist + 1;
		diagonal 	= dist + s2.length() - s1.length();

		matArr.resize(row);
		for(i = 0; i < row; ++i)
			matArr[i].resize(col, 0);

		//if(procID == 1 && checkTemp == 3164)
			//	cout << str1 << " -- " << str2 << " rt " << dist << endl;

		for(i = 0; i < dist + 1; i++)
		{
			for(j = dist - i; j < col; j++)
			{
				if (i == 0)
					matArr[i][j]	= j - dist;
				else if(j == (dist - i))
					matArr[i][j] 	= matArr[i - 1][j + 1] + 1;
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}

				if((j == diagonal) && matArr[i][j] > dist)
					return threshold + 1;
			}
		}

		for(i = dist + 1; i < s2.length() - dist + 1; i++)
		{
			for(j = 0; j < col; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else if(j != (col - 1))
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		for(i = s2.length() - dist + 1; i < row; i++)
		{
			for(j = 0; j < col - i + s2.length() - dist; j++)
			{
				if(j == 0)
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else
				{
					if((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return threshold + 1;
			}
		}

		//if(procID == 1 && checkTemp == 3164)
			//cout << str1 << " -- " << str2 << " hukjhk " << matArr[row - 1][diagonal] << endl;

		return matArr[row - 1][diagonal];

	}
}


int calculateBasicED2(string& str1, string& str2, int threshRem)
{
	int row, col, i, j;
	vector<vector<int> > matArr;

	row		= str1.length() + 1;
	col 	= str2.length() + 1;

	matArr.resize(row);
	for(i = 0; i < row; ++i)
		matArr[i].resize(col, 0);

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			if (i == 0)
				matArr[i][j] = j;
			else if (j == 0)
				matArr[i][j] = i;
			else
			{
				if((int)str1[i-1] == (int)str2[j-1])
					matArr[i][j]	= min(min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				else
					matArr[i][j] 	= min(min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
			}

			if((row - col) == (i - j) && (matArr[i][j] > threshRem))
			{
				return threshold + 1;
			}
		}
	}

	return (matArr[row-1][col-1]);
}



// change to lowercase
string convertToLower(string& str)
{
	for(unsigned int i = 0; i < str.length(); ++i)
		str[i]	= tolower(str[i]);
	return str;
}

// find root of a point in components

int findRoot(int pointID, vector<int> &parentArr)
{
	if (parentArr[pointID] != pointID)
		parentArr[pointID]	= findRoot(parentArr[pointID], parentArr);

	return parentArr[pointID];
}

// unify two components rooted at rootU and rootV

void makeEquivalent(int rootU, int rootV, vector<int> &parentArr, vector<int> &weightArr)
{
	if(weightArr[rootU] < weightArr[rootV])
		parentArr[rootU]	= rootV;
	else if(weightArr[rootU] > weightArr[rootV])
		parentArr[rootV]	= rootU;
	else
	{
		parentArr[rootV]	= rootU;
		weightArr[rootU]	= weightArr[rootU] + 1;
	}
}

