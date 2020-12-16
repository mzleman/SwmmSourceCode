//-----------------------------------------------------------------------------
//   input.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14  (Build 5.1.001)
//             09/15/14  (Build 5.1.007)
//             08/01/16  (Build 5.1.011)
//   Author:   L. Rossman
//
//   Input data processing functions.
//
//   Build 5.1.007:
//   - Support added for climate adjustment input data.
//
//   Build 5.1.011:
//   - Support added for reading hydraulic event dates.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"
#include "lid.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const int MAXERRS = 100;        // Max. input errors reported

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input
static int  Mobjects[MAX_OBJ_TYPES];   // Working number of objects of each type
static int  Mnodes[MAX_NODE_TYPES];    // Working number of node objects
static int  Mlinks[MAX_LINK_TYPES];    // Working number of link objects
static int  Mevents;                   // Working number of event periods

//-----------------------------------------------------------------------------
//  External Functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  input_countObjects  (called by swmm_open in swmm5.c)
//  input_readData      (called by swmm_open in swmm5.c)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int  addObject(int objType, char* id);
static int  getTokens(char *s);
static int  parseLine(int sect, char* line);
static int  readOption(char* line);
static int  readTitle(char* line);
static int  readControl(char* tok[], int ntoks);
static int  readNode(int type);
static int  readLink(int type);
static int  readEvent(char* tok[], int ntoks);

//=============================================================================

int input_countObjects()
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads input file to determine number of system objects.
//
{
    char  line[MAXLINE+1];             // line from input data file     
    char  wLine[MAXLINE+1];            // working copy of input line   
    char  *tok;                        // first string token of line          
    int   sect = -1, newsect;          // input data sections          
    int   errcode = 0;                 // error code
    int   errsum = 0;                  // number of errors found                   
    int   i;
    long  lineCount = 0;

    // --- initialize number of objects & set default values
    if ( ErrorCode ) return ErrorCode;
    error_setInpError(0, "");
    for (i = 0; i < MAX_OBJ_TYPES; i++) Nobjects[i] = 0;
    for (i = 0; i < MAX_NODE_TYPES; i++) Nnodes[i] = 0;
    for (i = 0; i < MAX_LINK_TYPES; i++) Nlinks[i] = 0;

    // --- make pass through data file counting number of each object
	/*如果文件中的该行，不足n-1个字符，则读完该行就结束。如若该行（包括最后一个换行符）的字符数超过n-1，则fgets只返回一个不完整的行，
	但是，缓冲区总是以NULL字符结尾，对fgets的下一次调用会继续读该行。
	函数调用失败或者读到文件尾返回null*/
    while ( fgets(line, MAXLINE, Finp.file) != NULL )
    {
        // --- skip blank lines & those beginning with a comment
        lineCount++;
        strcpy(wLine, line);   // make working copy of line
        tok = strtok(wLine, SEPSTR);   // get first text token on line
		//SEPSTR="\t\n\r",用这三个空字符取分割一行字符串，返回一组字符串 ，每次调用该函数 返回的数据类型是 char* ！！
		//tok指向切割后的第一个子串
		//第一次调用时返回切割后的第一个字符串
		//继续切割需要调用strtok(null, SEPSTR)，从上次切割位置继续


        if ( tok == NULL ) continue;//空行，则跳过对本行的读取
        if ( *tok == ';' ) continue;//开头为";",跳过对本行的读取

        // --- check if line begins with a new section heading
        if ( *tok == '[' )//读取到一个新的模块
        {
            // --- look for heading in list of section keywords
            newsect = findmatch(tok, SectWords);//返回匹配的在SectWords的索引，-1代表没找到，即发生错误
            if ( newsect >= 0 )//是一个合法的新模块名 
            {
                sect = newsect;//将sect改为newsect 
                continue;//跳过本循环，读取一下行文本，开始读取本模块的数据
            }
            else//模块名不合法
            {
                sect = -1;
                errcode = ERR_KEYWORD;
            }
        }

        // --- if in OPTIONS section then read the option setting
        //     otherwise add object and its ID name (tok) to project
		//一旦sect改变为某个整数，则接下来一直会按照该类型的数据读取每一非空行，直到sect改变为另外一个值
        if ( sect == s_OPTION ) errcode = readOption(line);//sect==1代表option 读取option，

        else if ( sect >= 0 )   errcode = addObject(sect, tok);//为sect模块添加id为tok的数据

        // --- report any error found
        if ( errcode )
        {
            report_writeInputErrorMsg(errcode, sect, line, lineCount);
            errsum++;
            if (errsum >= MAXERRS ) break;
        }
    }

    // --- set global error code if input errors were found
    if ( errsum > 0 ) ErrorCode = ERR_INPUT;
    return ErrorCode;
}

//=============================================================================

int input_readData()
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads input file to determine input parameters for each object.
//
{
    char  line[MAXLINE+1];        // line from input data file
    char  wLine[MAXLINE+1];       // working copy of input line
    char* comment;                // ptr. to start of comment in input line
    int   sect, newsect;          // data sections
    int   inperr, errsum;         // error code & total error count
    int   lineLength;             // number of characters in input line
    int   i;
    long  lineCount = 0;

    // --- initialize working item count arrays
    //     (final counts in Mobjects, Mnodes & Mlinks should
    //      match those in Nobjects, Nnodes and Nlinks).
    if ( ErrorCode ) return ErrorCode;

    error_setInpError(0, "");//初始化ErrString为""

    for (i = 0; i < MAX_OBJ_TYPES; i++)  Mobjects[i] = 0;
    for (i = 0; i < MAX_NODE_TYPES; i++) Mnodes[i] = 0;
    for (i = 0; i < MAX_LINK_TYPES; i++) Mlinks[i] = 0;
    Mevents = 0;

    // --- initialize starting date for all time series 初始化所有时间序列
    for ( i = 0; i < Nobjects[TSERIES]; i++ )
    {
        Tseries[i].lastDate = StartDate + StartTime;
    }

    // --- read each line from input file
    sect = 0;
    errsum = 0;
    rewind(Finp.file);//回到文件头部
    while ( fgets(line, MAXLINE, Finp.file) != NULL )
    {
        // --- make copy of line and scan for tokens
        lineCount++;
        strcpy(wLine, line);
        Ntokens = getTokens(wLine);

        // --- skip blank lines and comments
        if ( Ntokens == 0 ) continue;
        if ( *Tok[0] == ';' ) continue;

        // --- check if max. line length exceeded
        lineLength = strlen(line);
        if ( lineLength >= MAXLINE )
        {
            // --- don't count comment if present
            comment = strchr(line, ';'); //";"为一行中注释的标记
            if ( comment ) lineLength = comment - line;    // Pointer math here  ";"位置 - 行首位置 = ";"前的字符串长度
            if ( lineLength >= MAXLINE )
            {
				//在report文件中记录一个错误，行内字符串长度超过允许的最大长度
                inperr = ERR_LINE_LENGTH;
                report_writeInputErrorMsg(inperr, sect, line, lineCount);
                errsum++;
            }
        }

        // --- check if at start of a new input section
        if (*Tok[0] == '[')
        {
            // --- match token against list of section keywords
            newsect = findmatch(Tok[0], SectWords);
            
			if (newsect >= 0)//发现一个新模块
            {
                // --- SPECIAL CASE FOR TRANSECTS
                //     finish processing the last set of transect data

                if ( sect == s_TRANSECT ) //断面形状模块
                    transect_validate(Nobjects[TRANSECT]-1);

                // --- begin a new input section

                sect = newsect;
                continue;
            }
            else //模块名出错
            {
                inperr = error_setInpError(ERR_KEYWORD, Tok[0]);
                report_writeInputErrorMsg(inperr, sect, line, lineCount);
                errsum++;
                break;
            }
        }

        // --- otherwise parse tokens from input line
        else
        {
            inperr = parseLine(sect, line);
            if ( inperr > 0 )
            {
                errsum++;
                if ( errsum > MAXERRS ) report_writeLine(FMT19);
                else report_writeInputErrorMsg(inperr, sect, line, lineCount);
            }
        }

        // --- stop if reach end of file or max. error count
        if (errsum > MAXERRS) break;
    }   /* End of while */

    // --- check for errors
    if (errsum > 0)  ErrorCode = ERR_INPUT;
    return ErrorCode;
}

//=============================================================================

int  addObject(int objType, char* id)
//
//  Input:   objType = object type index
//           id = object's ID string
//  Output:  returns an error code
//  Purpose: adds a new object to the project.
//
{	/*----------------------------一共有16种Object--------------------------*/
    int errcode = 0;
    switch( objType )
    {
		/*
		int   project_addObject(int type, char *id, int n)会返回三种值

		 返回 1: 添加成功，在hash表中插入key和data，key为对象字符串名，data为对象在模块中是第几个
		 
		 返回-1:hash失败，但接下来不会报错，而且该模块的记录数+1 Nobjects[type]++;

		 返回 0：报错。
		*/

		// 3--代表雨量计模块
      case s_RAINGAGE:
        if ( !project_addObject(GAGE, id, Nobjects[GAGE]) )//0代表已经存在该id的数据 1代表在hash表中成功插入一条数据 -1代表hash失败或内存分配失败
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[GAGE]++;
        break;

		//  6--子汇水区模块
      case s_SUBCATCH:
        if ( !project_addObject(SUBCATCH, id, Nobjects[SUBCATCH]) )//内部会再次调用project_addObject(SUBCATCH, id, Nobjects[SUBCATCH])，如果发现会返回0 即出现error
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[SUBCATCH]++;//如果对id字符串hash失败，会出现异常 即该模块的数据个数算上了这一条内容，但hash表中没有插入这条记录
        break;

		//  9--含水层
      case s_AQUIFER:
        if ( !project_addObject(AQUIFER, id, Nobjects[AQUIFER]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[AQUIFER]++;
        break;

		// 34--单位过程线（单位量曲线）
      case s_UNITHYD:
        // --- the same Unit Hydrograph can span several lines//同一条曲线可能跨多行文本
        if ( project_findObject(UNITHYD, id) < 0 )
        {
            if ( !project_addObject(UNITHYD, id, Nobjects[UNITHYD]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            Nobjects[UNITHYD]++;
        }
        break;

		// 11--融雪
      case s_SNOWMELT:
        // --- the same Snowmelt object can appear on several lines//同一条曲线可能跨多行文本
        if ( project_findObject(SNOWMELT, id) < 0 )
        {
            if ( !project_addObject(SNOWMELT, id, Nobjects[SNOWMELT]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            Nobjects[SNOWMELT]++;
        }
        break;

		// 12--雨水井
      case s_JUNCTION: 
        if ( !project_addObject(NODE, id, Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[NODE]++; //节点个数+1
        Nnodes[JUNCTION]++; //雨水井个数+1
        break;

		// 13--出水口
      case s_OUTFALL:
        if ( !project_addObject(NODE, id, Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[NODE]++;//节点个数+1
        Nnodes[OUTFALL]++;//出水口个数+1
        break;

		// 14--蓄水设施
      case s_STORAGE:
        if ( !project_addObject(NODE, id, Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[NODE]++;//节点个数+1
        Nnodes[STORAGE]++;//蓄水设施节点个数+1
        break;

		// 15--分流器
      case s_DIVIDER: 
        if ( !project_addObject(NODE, id, Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[NODE]++;//节点个数+1
        Nnodes[DIVIDER]++;//分流器节点个数+1
        break;
		
		// 16--排水管
      case s_CONDUIT:
        if ( !project_addObject(LINK, id, Nobjects[LINK]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[LINK]++;//管线数+1
        Nlinks[CONDUIT]++;//排水管数+1
        break;

		// 17--水泵
      case s_PUMP:
        if ( !project_addObject(LINK, id, Nobjects[LINK]) ) 
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[LINK]++;//管线数+1
        Nlinks[PUMP]++;//水泵数+1
        break;
		
		// 18--孔(节流孔?)
      case s_ORIFICE:
        if ( !project_addObject(LINK, id, Nobjects[LINK]) ) 
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[LINK]++;//管线数+1
        Nlinks[ORIFICE]++;//节流空数+1
        break;

		// 19--堰
      case s_WEIR:
        if ( !project_addObject(LINK, id, Nobjects[LINK]) ) 
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[LINK]++;//管线数+1
        Nlinks[WEIR]++;//堰数+1
        break;

		// 20--排放孔
      case s_OUTLET:
        if ( !project_addObject(LINK, id, Nobjects[LINK]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[LINK]++;//管线数+1
        Nlinks[OUTLET]++;//排放孔数+1
        break;

		// 25--污染物
      case s_POLLUTANT:
        if ( !project_addObject(POLLUT, id, Nobjects[POLLUT]) ) 
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[POLLUT]++;//污染物数+1
        break;

		//  26--土地利用？？？
      case s_LANDUSE:
        if ( !project_addObject(LANDUSE, id, Nobjects[LANDUSE]) ) 
            errcode = error_setInpError(ERR_DUP_NAME, id);
        Nobjects[LANDUSE]++;
        break;

		// 32--时间模式？？
      case s_PATTERN:
        // --- a time pattern can span several lines
        if ( project_findObject(TIMEPATTERN, id) < 0 )
        {
            if ( !project_addObject(TIMEPATTERN, id, Nobjects[TIMEPATTERN]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            Nobjects[TIMEPATTERN]++;
        }
        break;
		
		// 37--曲线
      case s_CURVE:
        // --- a Curve can span several lines
        if ( project_findObject(CURVE, id) < 0 )
        {
            if ( !project_addObject(CURVE, id, Nobjects[CURVE]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            Nobjects[CURVE]++;

            // --- check for a conduit shape curve
            id = strtok(NULL, SEPSTR);
            if ( findmatch(id, CurveTypeWords) == SHAPE_CURVE )
                Nobjects[SHAPE]++;
        }
        break;

		// 38--时间序列
      case s_TIMESERIES:
        // --- a Time Series can span several lines
        if ( project_findObject(TSERIES, id) < 0 )
        {
            if ( !project_addObject(TSERIES, id, Nobjects[TSERIES]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            Nobjects[TSERIES]++;
        }
        break;

		// 24--控制
      case s_CONTROL:
        if ( match(id, w_RULE) ) Nobjects[CONTROL]++;
        break;

		// 22--横断面
      case s_TRANSECT:
        // --- for TRANSECTS, ID name appears as second entry on X1 line
        if ( match(id, "X1") )
        {
            id = strtok(NULL, SEPSTR);
            if ( id ) 
            {
                if ( !project_addObject(TRANSECT, id, Nobjects[TRANSECT]) )
                    errcode = error_setInpError(ERR_DUP_NAME, id);
                Nobjects[TRANSECT]++;
            }
        }
        break;

		// 49--低影响开发控制
      case s_LID_CONTROL:
        // --- an LID object can span several lines
        if ( project_findObject(LID, id) < 0 )
        {
            if ( !project_addObject(LID, id, Nobjects[LID]) )
            {
                errcode = error_setInpError(ERR_DUP_NAME, id);
            }
            Nobjects[LID]++;
        }
        break;

		// 53--事件
      case s_EVENT: 
		  NumEvents++; 
		  break;
    }
    return errcode;
}

//=============================================================================

int  parseLine(int sect, char *line)
//
//  Input:   sect  = current section of input file
//           *line = line of text read from input file
//  Output:  returns error code or 0 if no error found
//  Purpose: parses contents of a tokenized line of text read from input file.
//
{
    int j, err;
	int id;         //mz
	double x, y;	//mz
    switch (sect)
    {
      case s_TITLE:
        return readTitle(line);

      case s_RAINGAGE:
        j = Mobjects[GAGE];
        err = gage_readParams(j, Tok, Ntokens);
        Mobjects[GAGE]++;
        return err;

      case s_TEMP:
        return climate_readParams(Tok, Ntokens);

      case s_EVAP:
        return climate_readEvapParams(Tok, Ntokens);

      case s_ADJUST:
        return climate_readAdjustments(Tok, Ntokens);

      case s_SUBCATCH:
        j = Mobjects[SUBCATCH];
        err = subcatch_readParams(j, Tok, Ntokens);
        Mobjects[SUBCATCH]++;
        return err;

      case s_SUBAREA:
        return subcatch_readSubareaParams(Tok, Ntokens);

      case s_INFIL:
        return infil_readParams(InfilModel, Tok, Ntokens);

      case s_AQUIFER:
        j = Mobjects[AQUIFER];
        err = gwater_readAquiferParams(j, Tok, Ntokens);
        Mobjects[AQUIFER]++;
        return err;

      case s_GROUNDWATER:
        return gwater_readGroundwaterParams(Tok, Ntokens);

      case s_GWF:
        return gwater_readFlowExpression(Tok, Ntokens);

      case s_SNOWMELT:
        return snow_readMeltParams(Tok, Ntokens);

      case s_JUNCTION:
        return readNode(JUNCTION);

      case s_OUTFALL:
        return readNode(OUTFALL);

      case s_STORAGE:
        return readNode(STORAGE);

      case s_DIVIDER:
        return readNode(DIVIDER);

      case s_CONDUIT:
        return readLink(CONDUIT);

      case s_PUMP:
        return readLink(PUMP);

      case s_ORIFICE:
        return readLink(ORIFICE);

      case s_WEIR:
        return readLink(WEIR);

      case s_OUTLET:
        return readLink(OUTLET);

      case s_XSECTION:
        return link_readXsectParams(Tok, Ntokens);

      case s_TRANSECT:
        return transect_readParams(&Mobjects[TRANSECT], Tok, Ntokens);

      case s_LOSSES:
        return link_readLossParams(Tok, Ntokens);

      case s_POLLUTANT:
        j = Mobjects[POLLUT];
        err = landuse_readPollutParams(j, Tok, Ntokens);
        Mobjects[POLLUT]++;
        return err;

      case s_LANDUSE:
        j = Mobjects[LANDUSE];
        err = landuse_readParams(j, Tok, Ntokens);
        Mobjects[LANDUSE]++;
        return err;

      case s_BUILDUP:
        return landuse_readBuildupParams(Tok, Ntokens);

      case s_WASHOFF:
        return landuse_readWashoffParams(Tok, Ntokens);

      case s_COVERAGE:
        return subcatch_readLanduseParams(Tok, Ntokens);

      case s_INFLOW:
        return inflow_readExtInflow(Tok, Ntokens);

      case s_DWF:
        return inflow_readDwfInflow(Tok, Ntokens);

      case s_PATTERN:
        return inflow_readDwfPattern(Tok, Ntokens);

      case s_RDII:
        return rdii_readRdiiInflow(Tok, Ntokens);

      case s_UNITHYD:
        return rdii_readUnitHydParams(Tok, Ntokens);

      case s_LOADING:
        return subcatch_readInitBuildup(Tok, Ntokens);

      case s_TREATMENT:
        return treatmnt_readExpression(Tok, Ntokens);

      case s_CURVE:
        return table_readCurve(Tok, Ntokens);

      case s_TIMESERIES:
        return table_readTimeseries(Tok, Ntokens);

      case s_CONTROL:
        return readControl(Tok, Ntokens);

      case s_REPORT:
        return report_readOptions(Tok, Ntokens);

      case s_FILE:
        return iface_readFileParams(Tok, Ntokens);

      case s_LID_CONTROL:
        return lid_readProcParams(Tok, Ntokens);

      case s_LID_USAGE:
        return lid_readGroupParams(Tok, Ntokens);

      case s_EVENT:
        return readEvent(Tok, Ntokens);
	  case s_COORDINATE:
		id = project_findObject(NODE, Tok[0]);
		if (id < 0) return 0;
		if (!getDouble(Tok[1], &x))
			return error_setInpError(ERR_NUMBER, Tok[1]);
		if (!getDouble(Tok[2], &y))
			return error_setInpError(ERR_NUMBER, Tok[2]);
		Node[id].X = x;
		Node[id].Y = y;
		return 0;
      default: return 0;
    }
}

//=============================================================================

int readControl(char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads a line of input for a control rule.
//
{
    int index;
    int keyword;

    // --- check for minimum number of tokens
    if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");

    // --- get index of control rule keyword
    keyword = findmatch(tok[0], RuleKeyWords);
    if ( keyword < 0 ) return error_setInpError(ERR_KEYWORD, tok[0]);

    // --- if line begins a new control rule, add rule ID to the database
    if ( keyword == 0 )
    {
        if ( !project_addObject(CONTROL, tok[1], Mobjects[CONTROL]) )
        {
            return error_setInpError(ERR_DUP_NAME, Tok[1]);
        }
        Mobjects[CONTROL]++;
    }

    // --- get index of last control rule processed
    index = Mobjects[CONTROL] - 1;
    if ( index < 0 ) return error_setInpError(ERR_RULE, "");

    // --- add current line as a new clause to the control rule
    return controls_addRuleClause(index, keyword, Tok, Ntokens);
}

//=============================================================================

int readOption(char* line)
//
//  Input:   line = line of input data
//  Output:  returns error code
//  Purpose: reads an input line containing a project option.
//
{
    Ntokens = getTokens(line);//getTokens(line)应该是一个通用的读取每一行标识符、数据的功能函数
    if ( Ntokens < 2 ) return 0;//该行数据少于两个 即缺少option名或值
    return project_readOption(Tok[0], Tok[1]);//option模块，每行只包含key-value 
}

//=============================================================================

int readTitle(char* line)
//
//  Input:   line = line from input file
//  Output:  returns error code
//  Purpose: reads project title from line of input.
//
{
    int i, n;
    for (i = 0; i < MAXTITLE; i++)
    {
        // --- find next empty Title entry
        if ( strlen(Title[i]) == 0 )//在数组中找到一个空的位置
        {
            // --- strip line feed character from input line
            n = strlen(line);
            if (line[n-1] == 10) line[n-1] = ' ';//ASC码 10代表换行符 '\n',将换行符换成空格

            // --- copy input line into Title entry
            sstrncpy(Title[i], line, MAXMSG);//读取一行Title
            break;
        }
    }
    return 0;
}

//=============================================================================
    
int readNode(int type)
//
//  Input:   type = type of node
//  Output:  returns error code
//  Purpose: reads data for a node from a line of input.
//
{
    int j = Mobjects[NODE];
    int k = Mnodes[type];
    int err = node_readParams(j, type, k, Tok, Ntokens);//Tok保存了需要的内容 char* Tok[MAXTOKS]
    Mobjects[NODE]++;
    Mnodes[type]++;
    return err;
}

//=============================================================================

int readLink(int type)
//
//  Input:   type = type of link
//  Output:  returns error code
//  Purpose: reads data for a link from a line of input.
//
{
    int j = Mobjects[LINK];
    int k = Mlinks[type];
    int err = link_readParams(j, type, k, Tok, Ntokens);
    Mobjects[LINK]++;
    Mlinks[type]++;
    return err;
}

//=============================================================================

int  readEvent(char* tok[], int ntoks)
{
    DateTime x[4];

    if ( ntoks < 4 ) return error_setInpError(ERR_ITEMS, "");
    if ( !datetime_strToDate(tok[0], &x[0]) )
        return error_setInpError(ERR_DATETIME, tok[0]);
    if ( !datetime_strToTime(tok[1], &x[1]) )
        return error_setInpError(ERR_DATETIME, tok[1]);
    if ( !datetime_strToDate(tok[2], &x[2]) )
        return error_setInpError(ERR_DATETIME, tok[2]);
    if ( !datetime_strToTime(tok[3], &x[3]) )
        return error_setInpError(ERR_DATETIME, tok[3]);

    Event[Mevents].start = x[0] + x[1];
    Event[Mevents].end = x[2] + x[3];
    if ( Event[Mevents].start >= Event[Mevents].end )
       return error_setInpError(ERR_DATETIME, " - start date exceeds end date");
    Mevents++;
    return 0;
}

//=============================================================================

int  findmatch(char *s, char *keyword[])
//
//  Input:   s = character string
//           keyword = array of keyword strings
//  Output:  returns index of matching keyword or -1 if no match found  
//  Purpose: finds match between string and array of keyword strings.
//
{
   int i = 0;
   while (keyword[i] != NULL)
   {
      if (match(s, keyword[i])) return(i);
      i++;
   }
   return(-1);
}

//=============================================================================

int  match(char *str, char *substr)
//
//  Input:   str = character string being searched
//           substr = sub-string being searched for
//  Output:  returns 1 if sub-string found, 0 if not
//  Purpose: sees if a sub-string of characters appears in a string
//           (not case sensitive).
//
{
    int i,j;

    // --- fail if substring is empty
    if (!substr[0]) return(0);

    // --- skip leading blanks of str
    for (i = 0; str[i]; i++)
    {
        if (str[i] != ' ') break;
    }

    // --- check if substr matches remainder of str
    for (i = i,j = 0; substr[j]; i++,j++)
    {
        if (!str[i] || UCHAR(str[i]) != UCHAR(substr[j])) return(0);
    }
    return(1);
}

//=============================================================================

int  getInt(char *s, int *y)
//
//  Input:   s = a character string
//  Output:  y = converted value of s,
//           returns 1 if conversion successful, 0 if not
//  Purpose: converts a string to an integer number.
//
{
    double x;
    if ( getDouble(s, &x) )
    {
        if ( x < 0.0 ) x -= 0.01;
        else x += 0.01;
        *y = (int)x;
        return 1;
    }
    *y = 0;
    return 0;
}

//=============================================================================

int  getFloat(char *s, float *y)
//
//  Input:   s = a character string
//  Output:  y = converted value of s,
//           returns 1 if conversion successful, 0 if not
//  Purpose: converts a string to a single precision floating point number.
//
{
    char *endptr;
    *y = (float) strtod(s, &endptr);
    if (*endptr > 0) return(0);
    return(1);
}

//=============================================================================

int  getDouble(char *s, double *y)
//
//  Input:   s = a character string
//  Output:  y = converted value of s,
//           returns 1 if conversion successful, 0 if not
//  Purpose: converts a string to a double precision floating point number.
//
{
    char *endptr;
    *y = strtod(s, &endptr);//第二个参数类型为 char ** 用于修改某个指针的指向 该指针指向s转数字过程中不匹配的结束字符
    if (*endptr > 0) return(0);//如果结束的字符不是'\0',则有误
    return(1);
}

//=============================================================================

int  getTokens(char *s)
//
//  Input:   s = a character string
//  Output:  returns number of tokens found in s
//  Purpose: scans a string for tokens, saving pointers to them
//           in shared variable Tok[].
//
//  Notes:   Tokens can be separated by the characters listed in SEPSTR
//           (spaces, tabs, newline, carriage return) which is defined
//           in CONSTS.H. Text between quotes is treated as a single token.
//
{
    int  len, m, n;
    char *c;

    // --- begin with no tokens
    for (n = 0; n < MAXTOKS; n++) Tok[n] = NULL;
    n = 0;

    // --- truncate s at start of comment 
    c = strchr(s,';');//char *c指向s中第一个";"的位置
    if (c) *c = '\0';//如果一行中有";",则该位置为整行字符串的结尾；即分号是inp文件的注释符号
    len = strlen(s);//计算字符串长度，到第一个"\0"

    // --- scan s for tokens until nothing left
	//作用是
    while (len > 0 && n < MAXTOKS)
    {
		//C 库函数 size_t strcspn(const char *str1, const char *str2) 检索字符串 str1 开头连续有几个字符都不含字符串 str2 中的字符。
        m = strcspn(s,SEPSTR);              // find token length  计算当前循环下从s指向的当前位置开始 第一个标识符有多长

        if (m == 0) s++;      // no token found 如果第一个标识符长度为0 代表现在s指向的位置是空白字符

        else  //m>0代表找到了一个标识符字符串
        {
            if (*s == '"')                  // token begins with quote 如果第一个字符是引号 则去除引号的影响，提取引号内的内容
            {
                s++;                        // start token after quote 指针后移一位

                len--;          // reduce length of s  剩余长度-1，即减去一个引号的长度

                m = strcspn(s,"\"\n");      // find end quote or new line 以"和\n结尾，重新计算当前标识符长度（
											//效果是一遇到引号，后面的内容都认为是一个字符串的内容，直到右引号和换行符的出现）
            }
            s[m] = '\0';                    // null-terminate the token 将标识符后一个字符(空白字符或右引号)变为\0，
            Tok[n] = s;                     // save pointer to token 记录第n个标识符开始的位置，然后结尾为\0
            n++;                            // update token count
            s += m+1;                       // begin next token  即s指针移动到分隔符后一位
        }
        len -= m+1;                         // update length of s 剩余长度减去token长度m 再减去分隔符长度1
    }
    return(n);
}

//=============================================================================
