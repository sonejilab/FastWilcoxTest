/* source copied from https://github.com/stenver/wilcoxon-test*/

#include <Rcpp.h>

#include "WilcoxonTest.h"

WilcoxonTest::WilcoxonTest(float * _data, int _dataXsize, int _dataYsize, string _testIndexes, string _controlIndexes)
{
    data = _data;
    dataXsize = _dataXsize;
    dataYsize = _dataYsize;
    testIndexes = parseIntString(_testIndexes);
    controlIndexes = parseIntString(_controlIndexes);
    if(testIndexes->size() != controlIndexes->size()){
        cout << "Control and Test indexes are not the same size!" << endl;
        cout << "This program implements paired Wilcoxon test, so please use data groups with the same size." << endl;
        throw;
    }
    readApproximatePtable();
}

WilcoxonTest::WilcoxonTest(float * _data, int _dataXsize, int _dataYsize, vector<int> * _testIndexes, vector<int> * _controlIndexes)
{
    data = _data;
    dataXsize = _dataXsize;
    dataYsize = _dataYsize;
    testIndexes = _testIndexes;
    controlIndexes = _controlIndexes;
    if(testIndexes->size() != controlIndexes->size()){
        cout << "Control and Test indexes are not the same size!" << endl;
        cout << "This program implements paired Wilcoxon test, so please use data groups with the same size." << endl;
        throw;
    }
    readApproximatePtable();
}

vector<double> * WilcoxonTest::test()
{
    pValues = new vector<double>();
    int numberOfTests = testIndexes->size();
    for(int y = 0; y < numberOfTests; y++)
    { 
        vector<float> * absoluteValues = new vector<float>();
        vector<float> * signs = new vector<float>;
        float w = calculateWValue(y, absoluteValues, signs);
        int numberOfZeroes = getNumberOfZeroes(absoluteValues);
        pValues->push_back(calculatePValue(w, numberOfZeroes));
    }
    return pValues;
}

vector<int> * WilcoxonTest::parseIntString(string input)
{
    vector<int> * intVector = new vector<int>();
    vector<string> * intStrings = splitLine(input, ',');
    for(unsigned int i = 0; i < intStrings->size(); i++){
        intVector->push_back(atoi(intStrings->at(i).c_str()));
    }
    return intVector;
}

void WilcoxonTest::readApproximatePtable()
{
    string fileLocation = "../extdata/approximateTable500.txt";
    ifstream listFile(fileLocation.c_str());
    if (!listFile.is_open())
    {
        cerr << "P values table file does not exist at " << fileLocation << endl;
        throw;
    }
    approximatePTable = new std::vector<std::vector<approximatePosition> * >();
    
    string sLine = "";
    while (!listFile.eof())
    {
        getline(listFile, sLine);
        if(sLine.compare("") != 0)
        {
            approximatePTable->push_back(getPositions(sLine));
        }
    }
    listFile.close();
}

void WilcoxonTest::trim(string& str)
{
  string::size_type pos = str.find_last_not_of(' ');
  if(pos != string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != string::npos) str.erase(0, pos);
  }
  else str.erase(str.begin(), str.end());
}  

std::vector<approximatePosition> * WilcoxonTest::getPositions(string positionsLine)
{
    std::vector<approximatePosition> * positions = new std::vector<approximatePosition>();
    string trimmed_input = positionsLine.substr(1, positionsLine.length() - 1);
    std::vector<string> * positionsStrings = splitLine(trimmed_input, ']');
    for(unsigned int i = 0; i < positionsStrings->size(); i++){
        string rawPosistionString = positionsStrings->at(i);
        if(rawPosistionString.compare("") == 0) 
        {
          break;
        }
        int beginning = rawPosistionString.find('[');
        string warPositionStringFromBeginning = rawPosistionString.substr(beginning);
        int middle = warPositionStringFromBeginning.find(',');
        string x_str = warPositionStringFromBeginning.substr(1, middle - 1);
        string y_str = warPositionStringFromBeginning.substr(middle + 1, warPositionStringFromBeginning.length() - 1);
        trim(y_str);
        trim(x_str);
        approximatePosition pos;
        pos.x = atoi(x_str.c_str());
        pos.y = atof(y_str.c_str());
        positions->push_back(pos);
    }
    return positions;
}

std::vector<string> * WilcoxonTest::splitLine(string inputString, char lineSplit)
{
    std::istringstream stringStream ( inputString );
    std::vector<string> * datacells = new std::vector<string>();
    while (!stringStream.eof())
    {
        string newString;
        getline( stringStream, newString, lineSplit);
        datacells->push_back(newString);
    }
    return datacells;
}

float WilcoxonTest::calculateWValue(int yIndex, vector<float> * absoluteValues, vector<float> * signs)
{
    for (int i = 0; i < dataXsize; i++)
    {   
        int testIndex = testIndexes->at(yIndex);
        int controlIndex = controlIndexes->at(yIndex);
        float x1 = data[(testIndex * dataXsize) + i];
        float x2 = data[(controlIndex * dataXsize) + i];
        float value = x2-x1;
        absoluteValues->push_back(abs(value));
        signs->push_back(getSign(value));
    }

    //Sort the list
    quicksort(0, dataXsize-1, absoluteValues, signs);
     
    double * ranks = rankThePairs(yIndex, absoluteValues);
    float w = 0;
    for (int i = 0; i < dataXsize; i++)
    { 
        if (absoluteValues->at(i) != 0)
        {
           w += ranks[i] * signs->at(i);
        }
    }
    delete ranks;
    return abs(w);
}

double WilcoxonTest::calculatePValue(float w, int numberOfZeroes)
{
    int Nr = 0;
    for (int x = numberOfZeroes; x < dataXsize; x++)
    {
        Nr++;
    }
    float z = calculateZValue(w, Nr);
    if(Nr > 500){
        if(z < 0)
        {
            return gsl_cdf_gaussian_P(z, 1);
        }
        return gsl_cdf_gaussian_Q(z, 1);
    }
    else
    {
        return getApproximatePValue(w, z);
    }
}

float WilcoxonTest::calculateZValue(float w, int Nr)
{
    float sigma = sqrt((Nr * (Nr+1) * (2*Nr + 1)) / 6);
    return (w - 0.5) / sigma;
}

double WilcoxonTest::getApproximatePValue(float w, float z)
{ 
    std::vector<approximatePosition> * approximatePositions = approximatePTable->at(dataXsize);
    approximatePosition beginningPos = approximatePositions->at(0);
    for (unsigned int i = 1; i < approximatePositions->size(); i++)
    {   
        approximatePosition endPos = approximatePositions->at(i);
        if (w >= beginningPos.x && w <= endPos.x){
            return approximateP(w, z, beginningPos, endPos);
        }
        beginningPos = endPos;
    }
    return 0;
}

double WilcoxonTest::approximateP(float w, float z, approximatePosition beginningPos, approximatePosition endPos)
{
    double relativeValue; 
    if(w == beginningPos.x)
    {
        relativeValue = beginningPos.y;
    }
    else if (w == endPos.x)
    {
        relativeValue = endPos.y;
    }
    else
    {
        relativeValue = beginningPos.y + (endPos.y - beginningPos.y) * (abs(w) - beginningPos.x) / (endPos.x - beginningPos.x);
    }
    if (z < 0) {
      return relativeValue * gsl_cdf_gaussian_P(z, 1);
    }
    return relativeValue * gsl_cdf_gaussian_Q(z, 1); 
}

double * WilcoxonTest::rankThePairs(int yIndex, vector<float> * absoluteValues)
{
    double * ranks = new double[dataXsize];

    int i = 0;
    while (i < dataXsize)
    {
        int j = i + 1;
        while (j < dataXsize)
        {
            if(absoluteValues->at(i) != absoluteValues->at(j))
            {
                break;
            }
            j++;
        }
        for(int k = i; k <= j-1; k++)
        {   
            ranks[k] = 1 + (double)(i + j-1)/(double)2;
        }
        i = j;
    }
    return ranks;
}

int WilcoxonTest::getSign(float value)
{

    if(value > 0)
    {
        return 1;
    }
    else if(value < 0)
    {
        return -1;
    }
    return 0;
}

int WilcoxonTest::getNumberOfZeroes(vector<float> * absoluteValues)
{
    int numberOfZeroes = 0;
    for (int i = 0; i < dataXsize; i++)
    {
        if(absoluteValues->at(i) == 0)
        {
            numberOfZeroes++;
        }
    }
    return numberOfZeroes;
}

void WilcoxonTest::quicksort(int m, int n, vector<float> * absoluteValues, vector<float> * signs)
{
    float key;
    int i,j,k;
    if( m < n)
    {
        k = choose_pivot(m,n);
        swap(&absoluteValues->at(m), &absoluteValues->at(k));
        swap(&signs->at(m), &signs->at(k));
        key = absoluteValues->at(m);
        i = m + 1;
        j = n;
        while(i <= j)
        { 
            while((i <= n) && ( absoluteValues->at(i) <= key))
            {
                i++;
            }
            while((j >= m) && ( absoluteValues->at(j) > key))
            {
                j--;
            }
            if( i < j)
            {
        swap(&absoluteValues->at(i), &absoluteValues->at(j));
        swap(&signs->at(i), &signs->at(j));
            }
        }
        // swap two elements
        swap(&absoluteValues->at(m), &absoluteValues->at(j));
        swap(&signs->at(m), &signs->at(j));
        // recursively sort the lesser list
        quicksort(m,j - 1, absoluteValues, signs);
        quicksort(j + 1, n, absoluteValues, signs);
    }
}

void WilcoxonTest::swap(float *x,float *y)
{
    float temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

int WilcoxonTest::choose_pivot(int i,int j )
{
    return((i+j) /2);
}



