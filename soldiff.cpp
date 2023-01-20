#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

#ifndef PI
#define	PI 3.14159265358979
#endif

#ifndef D2R
#define D2R (PI/180.0)
#endif

#ifndef R2D
#define R2D (180.0/PI)
#endif

#ifndef TWOPI
#define	TWOPI (PI+PI)
#endif

#ifndef ae_WGS84
#define ae_WGS84 6378137.0
#endif

#ifndef finv_WGS84
#define finv_WGS84 298.257223563
#endif

#ifndef grav_WGS84
#define	grav_WGS84 9.7803267714e0
#endif

#define MAX_BUF_LEN (4096)

typedef struct
{
	uint8_t data[MAX_BUF_LEN];
	int nbyte;
	int nlen;
	int nmea_loc;
}parser_t;

double lat2local(double lat, double* lat2north)
{
	double f_WGS84 = (1.0 / finv_WGS84);
	double e2WGS84 = (2.0 * f_WGS84 - f_WGS84 * f_WGS84);
	double slat = sin(lat);
	double clat = cos(lat);
	double one_e2_slat2 = 1.0 - e2WGS84 * slat * slat;
	double Rn = ae_WGS84 / sqrt(one_e2_slat2);
	double Rm = Rn * (1.0 - e2WGS84) / (one_e2_slat2);
	*lat2north = Rm;
	return Rn * clat;
}

#define MAXFIELD 100

static int parse_fields(char* const buffer, char** val)
{
	char* p, *q;
	int n = 0;

	/* parse fields */
	for (p = buffer; *p && n < MAXFIELD; p = q + 1) {
		if (p == NULL) break;
		if ((q = strchr(p, ',')) || (q = strchr(p, '*')) || (q = strchr(p, '\n')) || (q = strchr(p, '\r'))) {
			val[n++] = p; *q = '\0';
		}
		else break;
	}
	return n;
}

static void set_output_file_name_without_directory(const char* fname, const char* key, char* outfname)
{
	char filename[255] = { 0 }, outfilename[255] = { 0 };
	strcpy(filename, fname);
	char* temp = strrchr(filename, '.');
	if (temp) temp[0] = '\0';
	temp = strrchr(filename, '/');
	if (temp)
	{
		sprintf(outfname, "%s%s", temp+1, key);
	}
	else
	{
		temp = strrchr(filename, '\\');
		if (temp)
		{
			sprintf(outfname, "%s%s", temp+1, key);
		}
		else
		{
			sprintf(outfname, "%s%s", filename, key);
		}
	}
}

static void set_output_file_name(const char* fname, const char* key, char *outfname)
{
	char filename[255] = { 0 }, outfilename[255] = { 0 };
	strcpy(filename, fname);
	char* temp = strrchr(filename, '.');
	if (temp) temp[0] = '\0';
	sprintf(outfname, "%s%s", filename, key);
}


static FILE* set_output_file(const char* fname, const char* key)
{
	char filename[255] = { 0 };
	set_output_file_name(fname, key, filename);
	return fopen(filename, "w");
}

static void deg2dms(double deg, double* dms)
{
	double sign = deg < 0.0 ? (-1.0) : (1.0), a = fabs(deg);
	dms[0] = floor(a); a = (a - dms[0]) * 60.0;
	dms[1] = floor(a); a = (a - dms[1]) * 60.0;
	dms[2] = a; dms[0] *= sign;
}
extern int outnmea_gga(unsigned char* buff, float time, int type, double* blh, int ns, float dop, float age)
{
	double h, ep[6], dms1[3], dms2[3];
	char* p = (char*)buff, * q, sum;

	if (type != 1 && type != 4 && type != 5) {
		p += sprintf(p, "$GPGGA,,,,,,,,,,,,,,");
		for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q;
		p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);
		return (int)(p - (char*)buff);
	}
	time -= 18.0;
	ep[2] = floor(time / (24 * 3600));
	time -= (float)(ep[2] * 24 * 3600.0);
	ep[3] = floor(time / 3600);
	time -= (float)(ep[3] * 3600);
	ep[4] = floor(time / 60);
	time -= (float)(ep[4] * 60);
	ep[5] = time;
	h = 0.0;
	deg2dms(fabs(blh[0]) * 180 / PI, dms1);
	deg2dms(fabs(blh[1]) * 180 / PI, dms2);
	p += sprintf(p, "$GPGGA,%02.0f%02.0f%06.3f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
		ep[3], ep[4], ep[5], dms1[0], dms1[1] + dms1[2] / 60.0, blh[0] >= 0 ? "N" : "S",
		dms2[0], dms2[1] + dms2[2] / 60.0, blh[1] >= 0 ? "E" : "W", type,
		ns, dop, blh[2] - h, h, age);
	for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q; /* check-sum */
	p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);
	return (int)(p - (char*)buff);
}

static int add_data_to_parser(parser_t* parser, int data)
{
	int ret = 0;
	if (parser->nbyte >= MAX_BUF_LEN) parser->nbyte = 0;
	if (parser->nbyte == 1 && data != 'G') parser->nbyte = 0;
	if (parser->nbyte == 0) memset(parser, 0, sizeof(parser_t));
	if (parser->nbyte == 0 && data != '$') return ret;
	parser->data[parser->nbyte++] = data;
	if (parser->nbyte > 3 && parser->data[parser->nbyte - 1] == 'G' && parser->data[parser->nbyte - 2] == '$')
	{
		memset(parser->data + 2, 0, sizeof(char) * (parser->nbyte - 2));
		parser->nbyte = 2;
	}
	if (parser->nbyte > 5 && data == '\r')
		ret = 0;
	if (parser->nbyte > 5 && parser->data[parser->nbyte - 1] == '\r' && parser->data[parser->nbyte - 4] == '*')
	{
		if (parser->nbyte<MAX_BUF_LEN)
			parser->data[parser->nbyte++] = '\n';
		//printf("%s", parser->data);
		ret = 1;
		parser->nlen = parser->nbyte;
		parser->nbyte = 0;
	}
	return ret;
}

typedef struct
{
	double time;
	double lat;
	double lon;
	double ht;
	double geod;
	double vn;
	double ve;
	double vu;
	double cdt;
	double cdt_drift;
	int nsat;
	int type;
	int age;
	int err;
	double cog;
	double pdop;
	double tdop;
	double hdop;
	double vdop;
	double hpl;
	double vpl;
	double speed;
	double heading;
	double acc_h;
	double acc_v;
	char ggatime[10];
}pvt_t;

static int read_nmea_gga(const char* fname, std::vector<pvt_t>& vpvt)
{
	FILE* fLOG = fopen(fname, "rb");
	if (!fLOG)
	{
		printf("cannot open %s\n", fname);
		return 0;
	}
	char buffer[512] = { 0 };
	FILE* fCSV = NULL;
	FILE* fGGA = NULL;
	FILE* fFIX = NULL;
	unsigned long index = 0;
	char* val[MAXFIELD];
	parser_t parser = { 0 };
	int data = 0;
	int ret = 0;
	while (fLOG != NULL && !feof(fLOG) && (data=fgetc(fLOG))!=EOF)
	{
		ret = add_data_to_parser(&parser, data);
		if (ret == 1)
		{
			strcpy(buffer, (char*)(parser.data));
			int num = parse_fields(buffer, val);
			if (num < 14 || strstr(val[0], "GGA") == NULL || strlen(val[3]) == 0 || strlen(val[5]) == 0) continue;
			pvt_t pvt = { 0 };
			strcpy(pvt.ggatime, val[1]);
			pvt.time = atof(val[1]);
			int hh = pvt.time / 10000;
			pvt.time = pvt.time - hh * 10000;
			int mm = pvt.time / 100;
			pvt.time = pvt.time - mm * 100;
			pvt.time = hh * 3600 + mm * 60 + pvt.time + 18.0;
			pvt.type = atoi(val[6]);
			pvt.nsat = atoi(val[7]);
			pvt.lat = atof(val[2]);
			int deg = pvt.lat / 100.0;
			pvt.lat = (pvt.lat - deg * 100);
			pvt.lat = (deg + pvt.lat / 60.0) * D2R;
			if (val[3][0] == 'S' || val[3][0] == 's')
				pvt.lat = -pvt.lat;
			pvt.lon = atof(val[4]);
			deg = pvt.lon / 100.0;
			pvt.lon = (pvt.lon - deg * 100);
			pvt.lon = (deg + pvt.lon / 60.0) * D2R;
			if (val[5][0] == 'W' || val[5][0] == 'w')
				pvt.lon = -pvt.lon;
			pvt.ht = atof(val[9]) + atof(val[11]);
			vpvt.push_back(pvt);
			if (!fGGA) fGGA = set_output_file(fname, "-gps.nmea");
			if (fGGA)
			{
				fprintf(fGGA, "%s", (char*)(parser.data));
			}
			if (pvt.type == 4)
			{
				if (!fFIX) fFIX = set_output_file(fname, "-fix.nmea");
				if (fFIX)
				{
					fprintf(fFIX, "%s", (char*)(parser.data));
				}
			}
			if (!fCSV) fCSV = set_output_file(fname, "-gps.csv");
			if (fCSV)
			{
				fprintf(fCSV, "%10.3f,%14.9f,%14.9f,%10.4f,%i,%i,%s\n", pvt.time, pvt.lat * R2D, pvt.lon * R2D, pvt.ht, pvt.nsat, pvt.type, pvt.ggatime);
			}
		}
	}
	if (fLOG) fclose(fLOG);
	if (fCSV) fclose(fCSV);
	if (fGGA) fclose(fGGA);
	if (fFIX) fclose(fFIX);
	return vpvt.size();
}

#define BAD_SAMPLE_INTERVAL 999999999

static int print_pvt_info(std::vector<pvt_t>& vpvt)
{
	unsigned long numoffix = 0;
	unsigned long numofrtk = 0;
	unsigned long numofspp = 0;
	unsigned long numofdgps = 0;
	unsigned long numofsol = 0;
	double avg[3] = { 0 };
	double avgfix[3] = { 0 };
	double avgrtk[3] = { 0 };
	double stdfix[3] = { 0 };
	double stdrtk[3] = { 0 };
	double std[3] = { 0 };
	std::vector<pvt_t>::iterator pvt = vpvt.begin();
	std::vector<pvt_t>::iterator pvt1= vpvt.begin();
	double l2e = 0;
	double l2n = 0;
	double dt = 0;
	double sample_interval = BAD_SAMPLE_INTERVAL;
	unsigned long numofsol_time_wrong = 0;
	unsigned long numofsol_datagap = 0;
	unsigned long numofepoch = 0;
	for (pvt = vpvt.begin(); pvt != vpvt.end(); ++pvt)
	{
		if (pvt->type == 1)
			++numofspp;
		else if (pvt->type == 2)
			++numofdgps;
		else if (pvt->type == 4)
		{
			++numoffix;
			++numofrtk;
			avgfix[0] += pvt->lat;
			avgfix[1] += pvt->lon;
			avgfix[2] += pvt->ht;
			avgrtk[0] += pvt->lat;
			avgrtk[1] += pvt->lon;
			avgrtk[2] += pvt->ht;
		}
		else if (pvt->type == 5)
		{
			++numofrtk;
			avgrtk[0] += pvt->lat;
			avgrtk[1] += pvt->lon;
			avgrtk[2] += pvt->ht;
		}
		avg[0] += pvt->lat;
		avg[1] += pvt->lon;
		avg[2] += pvt->ht;
	}
	/* fix solution stats */
	if (numoffix > 1)
	{
		avgfix[0] /= numoffix;
		avgfix[1] /= numoffix;
		avgfix[2] /= numoffix;
		for (pvt = vpvt.begin(); pvt != vpvt.end(); ++pvt)
		{
			if (pvt->type == 4)
			{
				stdfix[0] += (pvt->lat - avgfix[0]) * (pvt->lat - avgfix[0]);
				stdfix[1] += (pvt->lon - avgfix[1]) * (pvt->lon - avgfix[1]);
				stdfix[2] += (pvt->ht - avgfix[2]) * (pvt->ht - avgfix[2]);
			}
		}
		stdfix[0] = sqrt(stdfix[0] / (numoffix - 1));
		stdfix[1] = sqrt(stdfix[1] / (numoffix - 1));
		stdfix[2] = sqrt(stdfix[2] / (numoffix - 1));
		l2e = lat2local(avgfix[0], &l2n);
		stdfix[0] *= l2n;
		stdfix[1] *= l2e;
	}
	/* rtk solution stats */
	if (numofrtk > 1)
	{
		avgrtk[0] /= numofrtk;
		avgrtk[1] /= numofrtk;
		avgrtk[2] /= numofrtk;
		for (pvt = vpvt.begin(); pvt != vpvt.end(); ++pvt)
		{
			if (pvt->type == 4|| pvt->type == 5)
			{
				stdrtk[0] += (pvt->lat - avgrtk[0]) * (pvt->lat - avgrtk[0]);
				stdrtk[1] += (pvt->lon - avgrtk[1]) * (pvt->lon - avgrtk[1]);
				stdrtk[2] += (pvt->ht - avgrtk[2]) * (pvt->ht - avgrtk[2]);
			}
		}
		stdrtk[0] = sqrt(stdrtk[0] / (numofrtk - 1));
		stdrtk[1] = sqrt(stdrtk[1] / (numofrtk - 1));
		stdrtk[2] = sqrt(stdrtk[2] / (numofrtk - 1));
		l2e = lat2local(avgrtk[0], &l2n);
		stdrtk[0] *= l2n;
		stdrtk[1] *= l2e;
	}
	numofsol = vpvt.size();
	if (numofsol > 1)
	{
		avg[0] /= numofsol;
		avg[1] /= numofsol;
		avg[2] /= numofsol;
		for (pvt = vpvt.begin(); pvt != vpvt.end(); ++pvt)
		{
			std[0] += (pvt->lat - avg[0]) * (pvt->lat - avg[0]);
			std[1] += (pvt->lon - avg[1]) * (pvt->lon - avg[1]);
			std[2] += (pvt->ht - avg[2]) * (pvt->ht - avg[2]);
		}
		std[0] = sqrt(std[0] / (numofsol - 1));
		std[1] = sqrt(std[1] / (numofsol - 1));
		std[2] = sqrt(std[2] / (numofsol - 1));
		l2e = lat2local(avg[0], &l2n);
		std[0] *= l2n;
		std[1] *= l2e;
	}
	printf("number of sol [spp,dgps,rtk,fix,other] %i[%i,%i,%i,%i,%i], fix rate [%7.2f,%7.2f]\n", numofsol, numofspp, numofdgps, numofrtk, numoffix, numofsol - numofspp - numofdgps - numofrtk, (double)(numofsol > 0 ? (numoffix * 100.0) / numofsol : 0), (double)(numofrtk > 0 ? (numoffix * 100.0) / numofrtk : 0));
	if (numofsol > 1)
	{
		printf("sol [avg, std] %14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f\n", avg[0] * R2D, avg[1] * R2D, avg[2], std[0], std[1], std[2]);
	}
	if (numofrtk > 1)
	{
		printf("rtk [avg, std] %14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f\n", avgrtk[0] * R2D, avgrtk[1] * R2D, avgrtk[2], stdrtk[0], stdrtk[1], stdrtk[2]);
	}
	if (numoffix > 1)
	{
		printf("fix [avg, std] %14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f\n", avgfix[0] * R2D, avgfix[1] * R2D, avgfix[2], stdfix[0], stdfix[1], stdfix[2]);
	}
	if (numofrtk > 1)
	{
		l2e = lat2local(avgrtk[0], &l2n);
		printf("sol-rtk diff [m] %10.4f,%10.4f,%10.4f\n", (avg[0] - avgrtk[0])* l2n, (avg[1] - avgrtk[1])* l2e, (avg[2] - avgrtk[2]));
	}
	if (numoffix > 1)
	{
		l2e = lat2local(avgfix[0], &l2n);
		printf("rtk-fix diff [m] %10.4f,%10.4f,%10.4f\n", (avgrtk[0] - avgfix[0])* l2n, (avgrtk[1] - avgfix[1])* l2e, (avgrtk[2] - avgfix[2]));
	}
	/* found repated records and sample interval */
	for (pvt = vpvt.begin(); pvt != vpvt.end(); ++pvt)
	{
		if (pvt != vpvt.begin())
		{
			pvt1 = pvt - 1;
			dt = pvt->time - pvt1->time;
			if (dt < 0.0001)
			{
				++numofsol_time_wrong;
				printf("sol at time %10.4f has problems, with delt time %10.4f, total %i found with time problem\n", pvt1->time, dt, numofsol_time_wrong);
			}
			else
			{
				if (dt < sample_interval)
					sample_interval = dt;
			}
		}
	}
	if (sample_interval < BAD_SAMPLE_INTERVAL)
	{
		numofepoch = (vpvt.back().time - vpvt.front().time) / sample_interval + 1;
		printf("sample interval is %10.4f\n", sample_interval);
		printf("total epoch : %i\n", numofepoch);
		/* find data gaps */
		for (pvt = vpvt.begin(); pvt != vpvt.end(); ++pvt)
		{
			if (pvt != vpvt.begin())
			{
				pvt1 = pvt - 1;
				dt = pvt->time - pvt1->time;
				if (dt > sample_interval * (1.5))
				{
					++numofsol_datagap;
					printf("sol at time %10.4f has %10.4f data gap, total %i found with time problem %s\n", pvt1->time, dt, numofsol_datagap, pvt1->ggatime);
				}
			}
		}
		if (numofsol_datagap == 0)
		{
			printf("no data gap founded\n");
		}
		/* find RTK fix data gaps */
		pvt1 = vpvt.begin();
		unsigned long numoffix_cur = 0;
		for (pvt = vpvt.begin(); pvt != vpvt.end(); ++pvt)
		{
			if (pvt->type == 4)
			{
				/* RTK fix */
				if (numoffix_cur==0)
				{
					pvt1 = pvt;
				}
				else
				{
				}
				++numoffix_cur;
			}
			else
			{
				if (numoffix_cur > 0)
				{
					/* pre fix */
					dt = pvt->time - pvt1->time;
					printf("rtk fix at time %10.4f for %10.4f [s], GGA time at %s\n", pvt1->time, dt, pvt1->ggatime);
					numoffix_cur = 0;
				}
				else
				{

				}
			}
		}
		if (numoffix_cur > 0)
		{
			/* pre fix */
			dt = (vpvt.end()-1)->time - pvt1->time;
			printf("rtk fix at time %10.4f for %10.4f [s], GGA time at %s\n", pvt1->time, dt, pvt1->ggatime);
		}
	}

	return 0;
}

static int sol_diff(const char* fname1, const char* fname2, std::vector<pvt_t>& pvt1, std::vector<pvt_t>& pvt2)
{
	int numofmatch = 0;
	int i = 0;
	int j = 0;
	int wd1 = 0;
	int wd2 = 0;
	double time1 = 0;
	double time2 = 0;
	char fname1_buff[255] = { 0 };
	char fname2_buff[255] = { 0 };
	FILE* fDIF = NULL;
	set_output_file_name(fname1, "-", fname1_buff);
	set_output_file_name_without_directory(fname2, "-diff.csv", fname2_buff);
	std::vector<pvt_t>::iterator p1 = pvt1.begin();
	std::vector<pvt_t>::iterator p2 = pvt2.begin();

	for (; p1 != pvt1.end(); ++p1)
	{
		wd1 = p1->time / (24 * 3600);
		time1 = p1->time - wd1 * 24 * 3600;
		while (p2 != pvt2.end())
		{
			wd2 = p2->time / (24 * 3600);
			time2 = p2->time - wd2 * 24 * 3600;
			if (time1 > (time2 + 0.005))
			{
				++p2;
				continue;
			}
			break;
		}
		if (p2 == pvt2.end()) break;
		if (time1 < (time2 - 0.005))
			continue;
		double l2n = 0, l2e = lat2local(p1->lat, &l2n);
		double dn = (p1->lat - p2->lat) * l2n;
		double de = (p1->lon - p2->lon) * l2e;
		double dh = (p1->ht - p2->ht);
		if (!fDIF) fDIF = set_output_file(fname1_buff, fname2_buff);
		if (fDIF)
		{
			fprintf(fDIF, "%10.3f,%10.3f,%10.3f,%10.3f,%7.3f,%i,%i,%i,%i,%s\n", p1->time, dn, de, dh, time1 - time2, p1->type, p2->type, p1->nsat, p2->nsat, p1->ggatime);
		}
	}
	if (fDIF) fclose(fDIF);
	return numofmatch;
}
int main(int argc, char** argv)
{
	std::vector<pvt_t> vpvt1;
	std::vector<pvt_t> vpvt2;
	if (argc < 2)
	{
		/* */
		printf("%s filename1 filename2\n", argv[0]);
	}
	else
	{
		read_nmea_gga(argv[1], vpvt1);
		if (vpvt1.size() > 0)
		{
			printf("%s\n", argv[1]);
			print_pvt_info(vpvt1);
		}
		if (argc > 2)
		{
			read_nmea_gga(argv[2], vpvt2);
			if (vpvt2.size() > 0)
			{
				printf("%s\n", argv[2]);
				print_pvt_info(vpvt2);
			}
			sol_diff(argv[1], argv[2], vpvt1, vpvt2);
		}
	}
	return 0;
}