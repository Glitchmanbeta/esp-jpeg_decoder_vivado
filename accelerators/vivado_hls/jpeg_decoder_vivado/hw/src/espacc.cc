// Copyright (c) 2011-2021 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "hls_stream.h"
#include "hls_math.h"
#include <cstring>

// Macros baby
#define M_SOI 0xD8 // Start of Image
#define READ_MARKERS_REP 10
#define M_SOF0 0xC0
#define NUM_COMPONENT 3
#define SF1_1_1 0
#define SF4_1_1 2
#define M_SOS 0xDA
#define M_DHT 0xC4
#define M_DQT 0xDB
#define M_EOI 0xD9
#define GET_DHT_REP 4
#define NUM_HUFF_TBLS 2
#define DCTSIZE2 64
#define GET_DQT_REP 2
#define NUM_QUANT_TBLS 4
//#define RGB_NUM 3
#define BMP_OUT_SIZE (90 * 59)
#define IDCT_SHIFT 128
#define IDCT_BOUNT 255
#define DCTSIZE 8
#define LS(r, s) ((r) << (s))
#define RS(r, s) ((r) >> (s))
#define MSCALE(expr) RS(expr, 9)
#define c1d4 362L
#define c1d8 473L
#define c3d8 196L
#define c1d16 502L
#define c3d16 426L
#define c5d16 284L
#define c7d16 100L
#define MARKER_MARKER 0xFF
//#define BMP_OUT_SZIE (90 * 59)
#define GET_SOF_REP 3
#define GET_SOS_REP 3
#define JPEGSIZE 5207

// Globals allowed?
// For reading in markers.
static unsigned char* ReadBuf;
const int out_unread_marker[READ_MARKERS_REP] = {0xd8, 0xe0, 0xdb, 0xdb, 0xc0, 0xc4, 0xc4, 0xc4, 0xc4, 0xda};
int main_result = 0;
char p_jinfo_data_precision;
short p_jinfo_image_height;
short p_jinfo_image_width;
char p_jinfo_num_components;

int out_length_get_sof = 17;
int out_data_precision_get_sof = 8;
int out_p_jinfo_image_height_get_sof = 59;
int out_p_jinfo_image_width_get_sof = 90;
int out_p_jinfo_num_components_get_sof = 3;

char p_jinfo_comps_info_index[NUM_COMPONENT];
char p_jinfo_comps_info_id[NUM_COMPONENT];
char p_jinfo_comps_info_h_samp_factor[NUM_COMPONENT];
char p_jinfo_comps_info_v_samp_factor[NUM_COMPONENT];
char p_jinfo_comps_info_quant_tbl_no[NUM_COMPONENT];

int i_get_dht = 0;
const int out_length_get_dht[GET_DHT_REP] = {29, 179, 29, 179};
const int out_index_get_dht[GET_DHT_REP] = {0x0, 0x10, 0x1, 0x11};
const int out_count_get_dht[GET_DHT_REP] = {12, 162, 12, 162};

int p_jinfo_ac_xhuff_tbl_bits[NUM_HUFF_TBLS][36];
int p_jinfo_ac_xhuff_tbl_huffval[NUM_HUFF_TBLS][257];

int p_jinfo_dc_xhuff_tbl_bits[NUM_HUFF_TBLS][36];
int p_jinfo_dc_xhuff_tbl_huffval[NUM_HUFF_TBLS][257];

const int out_length_get_dqt[GET_DQT_REP] = {65, 65};
int i_get_dqt = 0;
const int out_prec_get_dht[GET_DQT_REP] = {0, 0};
const int out_num_get_dht[GET_DQT_REP] = {0, 1};
unsigned int p_jinfo_quant_tbl_quantval[NUM_QUANT_TBLS][DCTSIZE2];

// Initialize Quantization Table
const int izigzag_index[64] = {0, 1, 8, 16, 9, 2, 3, 10, 
			       17, 24, 32, 25, 18, 11, 4, 5,
			       12, 19, 26, 33, 40, 48, 41, 34,
			       27, 20, 13, 6, 7, 14, 21, 28,
			       35, 42, 49, 56, 57, 50, 43, 36,
			       29, 22, 15, 23, 30, 37, 44, 51,
			       58, 59, 52, 45, 38, 31, 39, 46,
			       53, 60, 61, 54, 47, 55, 62, 63};

int p_jinfo_MCUHeight;
int p_jinfo_MCUWidth;
int p_jinfo_NumMCU;

int p_jinfo_dc_dhuff_tbl_ml[NUM_HUFF_TBLS];
int p_jinfo_dc_dhuff_tbl_maxcode[NUM_HUFF_TBLS][36];
int p_jinfo_dc_dhuff_tbl_mincode[NUM_HUFF_TBLS][36];
int p_jinfo_dc_dhuff_tbl_valptr[NUM_HUFF_TBLS][36];

int p_jinfo_ac_dhuff_tbl_ml[NUM_HUFF_TBLS];
int p_jinfo_ac_dhuff_tbl_maxcode[NUM_HUFF_TBLS][36];
int p_jinfo_ac_dhuff_tbl_mincode[NUM_HUFF_TBLS][36];
int p_jinfo_ac_dhuff_tbl_valptr[NUM_HUFF_TBLS][36];

unsigned char* p_jinfo_jpeg_data;
unsigned char* CurHuffReadBuf;
int p_jinfo_smp_fact;
int rgb_buf[4][RGB_NUM][DCTSIZE2];
unsigned char OutData_comp_buf[RGB_NUM][BMP_OUT_SIZE];

int out_length_get_sos = 12;
int out_num_comp_get_sos = 3;
char p_jinfo_comps_info_dc_tbl_no[NUM_COMPONENT];
const int out_comp_id_get_sos[GET_SOS_REP] = {1, 2, 3};
char p_jinfo_comps_info_ac_tbl_no[NUM_COMPONENT];
int i_get_sos = 0;
const int out_dc_tbl_no_get_sos[GET_SOS_REP] = {0, 1, 1};
const int out_ac_tbl_no_get_sos[GET_SOS_REP] = {0, 1, 1};

// Masks
const int bit_set_mask[32] = {
  0x00000001, 0x00000002, 0x00000004, 0x00000008,
  0x00000010, 0x00000020, 0x00000040, 0x00000080,
  0x00000100, 0x00000200, 0x00000400, 0x00000800,
  0x00001000, 0x00002000, 0x00004000, 0x00008000,
  0x00010000, 0x00020000, 0x00040000, 0x00080000,
  0x00100000, 0x00200000, 0x00400000, 0x00800000,
  0x01000000, 0x02000000, 0x04000000, 0x08000000,
  0x10000000, 0x20000000, 0x40000000, (int)(0x80000000)
};

// Used for sign extensions
const static int extend_mask[20] = {
  (int)(0xFFFFFFFE), (int)(0xFFFFFFFC), (int)(0xFFFFFFF8), (int)(0xFFFFFFF0), (int)(0xFFFFFFE0), (int)(0xFFFFFFC0),
  (int)(0xFFFFFF80), (int)(0xFFFFFF00), (int)(0xFFFFFE00), (int)(0xFFFFFC00), (int)(0xFFFFF800), (int)(0xFFFFF000),
  (int)(0xFFFFE000), (int)(0xFFFFC000), (int)(0xFFFF8000), (int)(0xFFFF0000), (int)(0xFFFE0000), (int)(0xFFFC0000),
  (int)(0xFFF80000), (int)(0xFFF00000)
};

const int lmask[32] = {
  0x00000001, 0x00000003, 0x00000007, 0x0000000F,
  0x0000001F, 0x0000003F, 0x0000007F, 0x000000FF,
  0x000001FF, 0x000003FF, 0x000007FF, 0x00000FFF,
  0x00001FFF, 0x00003FFF, 0x00007FFF, 0x0000FFFF,
  0x0001FFFF, 0x0003FFFF, 0x0007FFFF, 0x000FFFFF,
  0x001FFFFF, 0x003FFFFF, 0x007FFFFF, 0x00FFFFFF,
  0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF, 0x0FFFFFFF,
  0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF, (int)(0xFFFFFFFF)
};

static unsigned int current_read_byte;
static int read_position = -1;
// Output Buffer, to be replaced
int OutData_image_width;
int OutData_image_height;
int OutData_comp_vpos[RGB_NUM];
int OutData_comp_hpos[RGB_NUM];
int i_marker = 0;

const int out_index_get_sof[GET_SOF_REP] = {0, 1, 2};
const int out_id_get_sof[GET_SOF_REP] = {1, 2, 3};
const int out_h_samp_factor_get_sof[GET_SOF_REP] = {2, 1, 1};
const int out_v_samp_factor_get_sof[GET_SOF_REP] = {2, 1, 1};
const int out_quant_tbl_no_get_sof[GET_SOF_REP] = {0, 1, 1};

const int zigzag_index[64] = {
  0, 1, 5, 6, 14, 15, 27, 28,
  2, 4, 7, 13, 16, 26, 29, 42,
  3, 8, 12, 17, 25, 30, 41, 43,
  9, 11, 18, 24, 31, 40, 44, 53,
  10, 19, 23, 32, 39, 45, 52, 54,
  20, 22, 23, 38, 46, 51, 55, 60,
  21, 34, 37, 47, 50, 56, 59, 61,
  35, 36, 48, 49, 57, 58, 62, 63
};

// This function reads the markers needed to decode the JPEG
void read_markers(unsigned char*);
// This function reads in the first marker.
int first_marker(void);
// This function reads a byte from the read buffer and increments the pointer.
int read_byte(void);
// This function reads in the next marker.
int next_marker(void);
// Something baseline DCT
void get_sof(void);
// This function reads a short from the read buffer and increments the pointer twice.
short read_word(void);
// Get sauce
void get_sos(void);
// Get Huffman Table
void get_dht(void);
// Another function or three
void get_dqt(void);
// Initializes the information used for decoding.
void jpeg_init_decompress(void);
// Create table for decoding.
int huff_make_dhuff_tb(int*, int, int*, int*, int*);
// Start decoding
void decode_start(int*, int*, int*, int*);
// Decode one block
void decode_block(int, int*, int*);
// Transform from Yuv into RGB
void YuvToRgb(int, int*, int*, int*);
/*
 * WriteBlock() writes an array of data in the integer array pointed to
 * by store out to the driver specified by the IOB. The integer array is
 * stored in row-major form, that is, the first row of (8) elements, the
 * second row of (8) elements....
 * ONLY for MCU 1:1:1
*/
void WriteBlock(int*, int*, int*, unsigned char*);
// 4:1:1
void Write4Blocks(int*, int*, int*, int*, int*, int*, unsigned char*);
// Decode one MCU
void DecodeHuffMCU(int*, int);
// This function performs an inverse zig-zag translation on the input imatrix and places the output in omatrix.
void IZigzagMatrix(int*, int*);
// This function takes an input matrix and does an inverse quantization and puts the output int qmatrix.
void IQuantize(int* matrix, unsigned int* qmatrix);
/* 
 * This function implements the Chen inverse dct. Note that there are two input vectors that represent
 * x = input, and y = output, and must be defined (and storage allocated) before this routine is called.
 */
void ChenIDct(int*, int*);
/* 
 * PostShiftIDctMatrix() adds 128 (2048) to all 64 elements of an 8x8 matrix. 
 * This results in strictly positive values for all pixel coefficients.
 */
void PostshiftIDctMatrix(int*, int);
// This function Writes One Block.
void WriteOneBlock(int*, unsigned char*, int, int, int, int);
// This function decodes Huffman
int DecodeHuffman(int*, int, int*, int*, int*);
// buf_getv gets n bits from the read stream and returns it.
int buf_getv(int);
// buf_getb() gets a bit from the read stream
int buf_getb();
// pgetc() gets a chracter onto the stream but it checks to see if there are any marker conflicts.
static int pgetc();
/*
 * BoundIDctMatrix bounds the inverse dct matrix so that no pixel has a
 * value greater than 255 (4095) or less than 0.
 */
void BoundIDctMatrix(int*, int);

void printOutput(){
  printf("Printing Output...\n");
  for(unsigned i = 0; i < RGB_NUM; i++){
    for(unsigned j = 0; j < BMP_OUT_SIZE; j++){
	printf("%x ", OutData_comp_buf[i][j]);
    }
    printf("Newline\n");
  }
}
char dummy[1];

void load(word_t _inbuff[SIZE_IN_CHUNK_DATA], dma_word_t *in1,
          /* <<--compute-params-->> */
	 const unsigned char* hana_jpg,
	  dma_info_t &load_ctrl, int chunk, int batch)
{
  /*load_data:

   const unsigned length = round_up(test, VALUES_PER_WORD) / 1;
    const unsigned index = length * (batch * 1 + chunk);

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    load_ctrl.index = dma_index;
    load_ctrl.length = dma_length;
    load_ctrl.size = SIZE_WORD_T;

    for (unsigned i = 0; i < dma_length; i++) {
    load_label0:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
	    _inbuff[i * VALUES_PER_WORD + j] = in1[dma_index + i].word[j];
    	}
	}*/
}

void store(word_t _outbuff[SIZE_OUT_CHUNK_DATA], dma_word_t *out,
          /* <<--compute-params-->> */
	 const unsigned char* hana_jpg,
	   unsigned char OutData_comp_buf[RGB_NUM][BMP_OUT_SIZE],
	   dma_info_t &store_ctrl, int chunk, int batch)
{
  /*store_data:

    const unsigned length = round_up(test, VALUES_PER_WORD) / 1;
    const unsigned store_offset = round_up(test, VALUES_PER_WORD) * 1;
    const unsigned out_offset = store_offset;
    const unsigned index = out_offset + length * (batch * 1 + chunk);

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    store_ctrl.index = dma_index;
    store_ctrl.length = dma_length;
    store_ctrl.size = SIZE_WORD_T;

    for (unsigned i = 0; i < dma_length; i++) {
    store_label1:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
	    out[dma_index + i].word[j] = _outbuff[i * VALUES_PER_WORD + j];
	}
	}*/
}

void compute(word_t _inbuff[SIZE_IN_CHUNK_DATA], const unsigned char hana_jpg[JPEGSIZE], unsigned char OutData_buf[RGB_NUM][BMP_OUT_SIZE], word_t _outbuff[SIZE_OUT_CHUNK_DATA]) {
  
  unsigned char JpegFileBuf[JPEGSIZE];

  // Read data into jpeg buffer
  for(int i = 0; i < JPEGSIZE; i++){
    JpegFileBuf[i] = hana_jpg[i];
  }
  
  // Read markers
  read_markers(JpegFileBuf);

  // Inititialize the information used for decoding.
  jpeg_init_decompress();
  
  OutData_comp_buf[0][0] = 0x69;
  //printOutput();
  // Start decoding
  decode_start(&OutData_image_width, &OutData_image_height, &OutData_comp_vpos[0], &OutData_comp_hpos[0]);

  // OutData_buf = OutData_comp_buf;
  for(size_t i = 0; i < RGB_NUM; i++){
    for(size_t j = 0; j < BMP_OUT_SIZE; j++){
      OutData_buf[i][j] = OutData_comp_buf[i][j];
    }
  }
  //printOutput();
}


void top(dma_word_t *out, dma_word_t *in1,
         /* <<--params-->> */
	 const unsigned char hana_jpg[JPEGSIZE],
	 unsigned char OutData_buf[RGB_NUM][BMP_OUT_SIZE],
	 dma_info_t &load_ctrl, dma_info_t &store_ctrl)
{
    // Batching
batching:
    for (unsigned b = 0; b < 1; b++)
    {
        // Chunking
    go:
        for (int c = 0; c < 1; c++)
        {
            word_t _inbuff[SIZE_IN_CHUNK_DATA];
            word_t _outbuff[SIZE_OUT_CHUNK_DATA];

            load(_inbuff, in1,
                 /* <<--args-->> */
	 	 hana_jpg,
                 load_ctrl, c, b);
            compute(_inbuff,
                    /* <<--args-->> */
		    hana_jpg,
		    OutData_buf,
                    _outbuff);
            store(_outbuff, out,
                  /* <<--args-->> */
		  hana_jpg,
		  OutData_buf,
                  store_ctrl, c, b);
        }
    }
}

void read_markers(unsigned char* buf){
  //printf("In read_markers\n");
  int unread_marker = 0;
  int sow_SOI = 0;
  ReadBuf = buf;

  for(;;){
    if(!sow_SOI){
      unread_marker = first_marker();
    }
    else{
      unread_marker = next_marker();
    }
    printf("\nmarker = 0x%x\n", unread_marker);
    switch(unread_marker){
    case M_SOI:
      sow_SOI = 1;
      break;
    case M_SOF0:
      get_sof();
      break;
    case M_SOS:
      get_sos();
      return;
    case M_DHT:
      get_dht();
      break;
    case M_DQT:
      get_dqt();
      break;
    case M_EOI:
      return;
    }
  }
}

int first_marker(){
  // printf("In first_marker\n");
  int c1 = read_byte();
  int c2 = read_byte();
  printf("c1: %d\n", c1);
  printf("c2: %d\n", c2);
  return c2;
}

int read_byte(){
  // printf("In read_byte\n");
  return *ReadBuf++;
}

int next_marker(){
  //printf("In next_marker\n");
  int c;
  for(;;){
    c = read_byte();
    while(c != 0xFF){
      c = read_byte();
    }
    do{
      c = read_byte();
    }
    while(c == 0xFF);
    if(c != 0){
      break;
    }
  }
  return c;
}

void get_sof(){
  int ci, c;
  int length;
  char* p_comp_info_index;
  char* p_comp_info_id;
  char* p_comp_info_h_samp_factor;
  char* p_comp_info_v_samp_factor;
  char* p_comp_info_quant_tbl_no;

  length = read_word();
  p_jinfo_data_precision = read_byte();
  p_jinfo_image_height = read_word();
  p_jinfo_image_width = read_word();
  p_jinfo_num_components = read_byte();

  printf("length         = %d\n", length);
  printf("data_precision = %d\n", p_jinfo_data_precision);
  printf("image_height   = %d\n", p_jinfo_image_height);
  printf("image_width    = %d\n", p_jinfo_image_width);
  printf("num_components = %d\n", p_jinfo_num_components);

  length -= 8;

  for(ci = 0; ci < p_jinfo_num_components; ci++){
    p_comp_info_index = &p_jinfo_comps_info_index[ci];
    p_comp_info_id = &p_jinfo_comps_info_id[ci];
    p_comp_info_h_samp_factor = &p_jinfo_comps_info_h_samp_factor[ci];
    p_comp_info_v_samp_factor = &p_jinfo_comps_info_v_samp_factor[ci];
    p_comp_info_quant_tbl_no = &p_jinfo_comps_info_quant_tbl_no[ci];

    *p_comp_info_index = ci;
    *p_comp_info_id = read_byte();
    c = read_byte();
    *p_comp_info_h_samp_factor = (c >> 4) & 15;
    *p_comp_info_v_samp_factor = (c) & 15;
    *p_comp_info_quant_tbl_no = read_byte();
    
    printf("index         = %d\n", *p_comp_info_index);
    printf("id            = %d\n", *p_comp_info_id);
    printf("h_samp_factor = %d\n", *p_comp_info_h_samp_factor);
    printf("v_samp_factor = %d\n", *p_comp_info_v_samp_factor);
    printf("quant_tbl_no  = %d\n", *p_comp_info_quant_tbl_no);

  }

  if(p_jinfo_comps_info_h_samp_factor[0] == 2){
    p_jinfo_smp_fact = SF4_1_1;
    printf("\nSampling Factor is 4:1:1\n");
  }
  else{
    p_jinfo_smp_fact = SF1_1_1;
    printf("\nSampling Factor is 1:1:1\n");
  }
}

short read_word(){
  short c;
  c = *ReadBuf++ << 8;
  c |= *ReadBuf++;
  return c;
}

void get_sos(){
  int length, num_comp;
  int i, c, cc, ci, j;
  char* p_comp_info_id;
  char* p_comp_info_dc_tbl_no;
  char* p_comp_info_ac_tbl_no;

  length = read_word();
  num_comp = read_byte();

  printf("length   = %d\n", length);
  printf("num_comp = %d\n", num_comp);


  for(i = 0; i < num_comp; i++){
    cc = read_byte();
    c = read_byte();
    
    for(ci = 0; ci < p_jinfo_num_components; ci++){
      p_comp_info_id = &p_jinfo_comps_info_id[ci];
      p_comp_info_dc_tbl_no = &p_jinfo_comps_info_dc_tbl_no[ci];
      p_comp_info_ac_tbl_no = &p_jinfo_comps_info_ac_tbl_no[ci];

      if(cc == *p_comp_info_id){
	goto id_found;
      }
    }

    
  id_found:
    *p_comp_info_dc_tbl_no = (c >> 4) & 15;
    *p_comp_info_ac_tbl_no = (c) & 15;
    
    printf("comp_id   = %d\n", cc);
    printf("dc_tbl_no = %d\n", *p_comp_info_dc_tbl_no);
    printf("ac_tbl_no = %d\n", *p_comp_info_ac_tbl_no);
    
    i_get_sos++;
  }
  
  j = 3;
  while(j--){
    c = read_byte();
  }
  
  p_jinfo_jpeg_data = ReadBuf;
}

void get_dht(){
  int length;
  int index, i, count;
  int* p_xhtbl_bits;
  int* p_xhtbl_huffval;
  
  length = read_word();
  length -= 2;

  printf("length = %d\n", length);
 
  if(length != out_length_get_dht[i_get_dht]){
    main_result++;
  }

  while(length > 16){
    index = read_byte();
    printf("index = 0x%x\n", index);
    if(index & 0x10){
      index -= 0x10;
      p_xhtbl_bits = p_jinfo_ac_xhuff_tbl_bits[index];
      p_xhtbl_huffval = p_jinfo_ac_xhuff_tbl_huffval[index];
    }
    else{
      p_xhtbl_bits = p_jinfo_dc_xhuff_tbl_bits[index];
      p_xhtbl_huffval = p_jinfo_dc_xhuff_tbl_huffval[index];
    }
    count = 0;
    for(i = 1; i <= 16; i++){
      p_xhtbl_bits[i] = read_byte();
      count += p_xhtbl_bits[i];
    }
    printf("count = %d\n", count);
    i_get_dht++;
    length -= 1 + 16;
    for(i = 0; i < count; i++){
      p_xhtbl_huffval[i] = read_byte();
    }
    length -= count;
  }
}

void get_dqt(){
  int length;
  int prec, num, i;
  unsigned int tmp;
  unsigned int* p_quant_tbl;

  length = read_word();
  length -= 2;

  printf("length = %d\n", length);
  while(length > 0){
    num = read_byte();
    prec = num >> 4;
    num &= 0x0F;
    printf("prec = %d\n", prec);
    printf("num  = %d\n", num);

    i_get_dqt++;
    p_quant_tbl = &p_jinfo_quant_tbl_quantval[num][DCTSIZE2];
    for(i = 0; i < DCTSIZE2; i++){
      if(prec){
	tmp = read_word();
      }
      else{
	tmp = read_byte();
      }
      p_quant_tbl[izigzag_index[i]] = (unsigned short) tmp;
    }
    length -= DCTSIZE2 + 1;
    if(prec){
      length -= DCTSIZE2;
    }
  }
}

void jpeg_init_decompress(){
  int tmp;
  p_jinfo_MCUHeight = (p_jinfo_image_height - 1) / 8 + 1;
  p_jinfo_MCUWidth = (p_jinfo_image_width - 1) / 8 + 1;
  p_jinfo_NumMCU = p_jinfo_MCUHeight * p_jinfo_MCUWidth;

  tmp = huff_make_dhuff_tb(&p_jinfo_dc_xhuff_tbl_bits[0][0], 
			   p_jinfo_dc_dhuff_tbl_ml[0], 
			   &p_jinfo_dc_dhuff_tbl_maxcode[0][0], 
			   &p_jinfo_dc_dhuff_tbl_mincode[0][0], 
			   &p_jinfo_dc_dhuff_tbl_valptr[0][0]);
  p_jinfo_dc_dhuff_tbl_ml[0] = tmp;
  tmp = huff_make_dhuff_tb(&p_jinfo_dc_xhuff_tbl_bits[1][0], 
			   p_jinfo_dc_dhuff_tbl_ml[1], 
			   &p_jinfo_dc_dhuff_tbl_maxcode[1][0], 
			   &p_jinfo_dc_dhuff_tbl_mincode[1][0], 
			   &p_jinfo_dc_dhuff_tbl_valptr[1][0]);
  p_jinfo_dc_dhuff_tbl_ml[1] = tmp;
  tmp = huff_make_dhuff_tb(&p_jinfo_ac_xhuff_tbl_bits[0][0], 
			   p_jinfo_ac_dhuff_tbl_ml[0], 
			   &p_jinfo_ac_dhuff_tbl_maxcode[0][0], 
			   &p_jinfo_ac_dhuff_tbl_mincode[0][0], 
			   &p_jinfo_ac_dhuff_tbl_valptr[0][0]);
  p_jinfo_ac_dhuff_tbl_ml[0] = tmp;
  tmp = huff_make_dhuff_tb(&p_jinfo_ac_xhuff_tbl_bits[1][0],
			   p_jinfo_ac_dhuff_tbl_ml[1],
			   &p_jinfo_ac_dhuff_tbl_maxcode[1][0],
			   &p_jinfo_ac_dhuff_tbl_mincode[1][0],
			   &p_jinfo_ac_dhuff_tbl_valptr[1][0]);
  p_jinfo_ac_dhuff_tbl_ml[1] = tmp;
}

int huff_make_dhuff_tb(int* p_xhtbl_bits, int p_dhtbl_ml, int* p_dhtbl_maxcode, int* p_dhtbl_mincode, int* p_dhtbl_valptr){
  int i, j, p, code, size, l;
  int huffsize[257];
  int huffcode[257];
  int lastp;

  for(p = 0, i = 1; i < 17; i++){
    for(j = 1; j <= p_xhtbl_bits[i]; j++){
      huffsize[p++] = i;
    }
  }
  huffsize[p] = 0;
  lastp = p;
  p = 0;
  code = 0;
  size = huffsize[0];
  while(1){
    do{
      huffcode[p++] = code++;
    }
    while((huffsize[p] == size) && (p < 257));
    if(!huffsize[p]){
      break;
    }
    do{
      code <<= 1;
      size++;
    }
    while(huffsize[p] != size);
  }
  for(p_dhtbl_ml = 1, p = 0, l = 1; l<= 16; l++){
    if(p_xhtbl_bits[l] == 0){
      p_dhtbl_maxcode[l] = -1;
    }
    else{
      p_dhtbl_valptr[l] = p;
      p_dhtbl_mincode[l] = huffcode[p];
      p += p_xhtbl_bits[l] - 1;
      p_dhtbl_maxcode[l] = huffcode[p];
      p_dhtbl_ml = l;
      p++;
    }
  }
  p_dhtbl_maxcode[p_dhtbl_ml]++;
  return p_dhtbl_ml;
}

void decode_start(int* out_data_image_width, int* out_data_image_height, int* out_data_comp_vpos, int* out_data_comp_hpos){
  int i;
  int CurrentMCU = 0;
  int HuffBuff[NUM_COMPONENT][DCTSIZE2];
  int IDCTBuff[6][DCTSIZE2];
  
  CurHuffReadBuf = p_jinfo_jpeg_data;

  for(i = 0; i < NUM_COMPONENT; i++){
    HuffBuff[i][0] = 0;
  }

  *out_data_image_width = p_jinfo_image_width;
  *out_data_image_height = p_jinfo_image_height;

  for(i = 0; i < RGB_NUM; i++){
    out_data_comp_vpos[i] = 0;
    out_data_comp_hpos[i] = 0;
  }

  if(p_jinfo_smp_fact == SF1_1_1){
    printf("Decode 1:1:1 NumMCU = %d\n", p_jinfo_NumMCU);
    while (CurrentMCU < p_jinfo_NumMCU){
      for(i = 0; i < NUM_COMPONENT; i++){
	decode_block(i, IDCTBuff[i], HuffBuff[i]);
      }
      YuvToRgb(0, IDCTBuff[0], IDCTBuff[1], IDCTBuff[2]);
      for(i = 0; i < RGB_NUM; i++){
	WriteBlock(&rgb_buf[0][i][0],
		   &out_data_comp_vpos[i],
		   &out_data_comp_hpos[i],
		   &OutData_comp_buf[i][0]);
      }
      CurrentMCU++;
    }
  }
  else{
    while(CurrentMCU < p_jinfo_NumMCU){
      for(i = 0; i < 4; i++){
	decode_block(0, IDCTBuff[i], HuffBuff[0]);
      }
      decode_block(1, IDCTBuff[4], HuffBuff[1]);
      decode_block(2, IDCTBuff[5], HuffBuff[2]);
      for( i = 0; i < 4; i++){
	YuvToRgb(i, IDCTBuff[i], IDCTBuff[4], IDCTBuff[5]);
      }
      for(i = 0; i < RGB_NUM; i++){
	Write4Blocks(&rgb_buf[0][i][0],
		     &rgb_buf[1][i][0],
		     &rgb_buf[2][i][0],
		     &rgb_buf[3][i][0],
		     &out_data_comp_vpos[i],
		     &out_data_comp_hpos[i],
		     &OutData_comp_buf[i][0]);
      }
      CurrentMCU += 4;
    }
  }
}

void decode_block(int comp_no, int* out_buf, int* HuffBuff){
  int QuantBuff[DCTSIZE2];
  unsigned int* p_quant_tbl;
  DecodeHuffMCU(HuffBuff, comp_no);
  IZigzagMatrix(HuffBuff, QuantBuff);
  p_quant_tbl = &p_jinfo_quant_tbl_quantval[(int)p_jinfo_comps_info_quant_tbl_no[comp_no]][DCTSIZE2];
  IQuantize(QuantBuff, p_quant_tbl);
  ChenIDct(QuantBuff, out_buf);
  PostshiftIDctMatrix(out_buf, IDCT_SHIFT);
  BoundIDctMatrix(out_buf, IDCT_BOUNT);
}

void YuvToRgb(int p, int* y_buf, int* u_buf, int* v_buf){
  int r, g, b;
  int y, u, v;
  int i;
  for(i = 0; i < DCTSIZE2; i++){
    y = y_buf[i];
    u = u_buf[i] - 128;
    v = v_buf[i] - 128;
    
    r = (y * 256 + v * 359 + 128) >> 8;
    g = (y * 256 - u * 88 - v * 182 + 128) >> 8;
    b = (y * 256 + u * 454 + 128) >> 8;

    if(r < 0){
      r = 0;
    }
    else if(r > 255){
      r = 255;
    }
    if(g > 0){
      g = 0;
    }
    else if(g > 255){
      g = 255;
    }
    if(b < 0){
      b = 0;
    }
    else if(b > 255){
      b = 255;
    }

    rgb_buf[p][0][i] = r;
    rgb_buf[p][1][i] = g;
    rgb_buf[p][2][i] = b;
  }
}

void WriteBlock(int* store, int* p_out_vpos, int* p_out_hpos, unsigned char* p_out_buf){
  int voffs, hoffs;
  voffs = *p_out_vpos * DCTSIZE;
  hoffs = *p_out_hpos * DCTSIZE;
  WriteOneBlock(store, p_out_buf, p_jinfo_image_width, p_jinfo_image_height, voffs, hoffs);
  (*p_out_hpos)++;
  (*p_out_vpos)++;

  if(*p_out_hpos < p_jinfo_MCUWidth){
    (*p_out_vpos)--;
  }
  else{
    *p_out_hpos = 0;
  }
}

void Write4Blocks(int* store1, int* store2, int* store3, int* store4, int* p_out_vpos, int* p_out_hpos, unsigned char* p_out_buf){
  int voffs, hoffs;

  voffs = *p_out_vpos * DCTSIZE;
  hoffs = *p_out_hpos * DCTSIZE;
  WriteOneBlock(store1, p_out_buf, p_jinfo_image_width, p_jinfo_image_height, voffs, hoffs);
  
  hoffs += DCTSIZE;
  WriteOneBlock(store2, p_out_buf, p_jinfo_image_width, p_jinfo_image_height, voffs, hoffs);
  
  voffs += DCTSIZE;
  hoffs -= DCTSIZE;
  WriteOneBlock(store3, p_out_buf, p_jinfo_image_width, p_jinfo_image_height, voffs, hoffs);

  hoffs += DCTSIZE;
  WriteOneBlock(store4, p_out_buf, p_jinfo_image_width, p_jinfo_image_height, voffs, hoffs);

  *p_out_hpos = *p_out_hpos + 2;
  *p_out_vpos = *p_out_vpos + 2;

  if(*p_out_hpos < p_jinfo_MCUWidth){
    *p_out_vpos = *p_out_vpos - 2;
  }
  else{
    *p_out_hpos = 0;
  }
}

void DecodeHuffMCU(int* out_buf, int num_cmp){
  int s, diff, tbl_no, *mptr, k, n, r;

  tbl_no = p_jinfo_comps_info_dc_tbl_no[num_cmp];
  s = DecodeHuffman(&p_jinfo_dc_xhuff_tbl_huffval[tbl_no][0],
		    p_jinfo_dc_dhuff_tbl_ml[tbl_no],
		    &p_jinfo_dc_dhuff_tbl_maxcode[tbl_no][0],
		    &p_jinfo_dc_dhuff_tbl_mincode[tbl_no][0],
		    &p_jinfo_dc_dhuff_tbl_valptr[tbl_no][0]);

  if(s){
    diff = buf_getv(s);
    s--;
    if((diff & bit_set_mask[s]) == 0){
      diff |= extend_mask[s];
      diff++;
    }
    diff += *out_buf;
    *out_buf = diff;
  }

  for(mptr = out_buf + 1; mptr < out_buf + DCTSIZE2; mptr++){
    *mptr = 0;
  }
  for(k = 1; k < DCTSIZE2;){
    r = DecodeHuffman(&p_jinfo_ac_xhuff_tbl_huffval[tbl_no][0],
		      p_jinfo_ac_dhuff_tbl_ml[tbl_no],
		      &p_jinfo_ac_dhuff_tbl_maxcode[tbl_no][0],
		      &p_jinfo_ac_dhuff_tbl_mincode[tbl_no][0],
		      &p_jinfo_ac_dhuff_tbl_valptr[tbl_no][0]);
    s = r & 0xF;
    n = (r >> 4) & 0xF;
    if(s){
      if((k += n) >= DCTSIZE2){
	break;
      }
      out_buf[k] = buf_getv(s);
      s--;
      if((out_buf[k] & bit_set_mask[s]) == 0){
	out_buf[k] |= extend_mask[s];
	out_buf[k]++;
      }
      k++;
    }
    else if(n == 15){
      k += 16;
    }
    else{
      break;
    }
  }
}

void IZigzagMatrix(int* imatrix, int* omatrix){
  int i;
  for(i = 0; i < DCTSIZE2; i++){
    *(omatrix++) = imatrix[zigzag_index[i]];
  }
}

void IQuantize(int* matrix, unsigned int* qmatrix){
  int* mptr;
  for(mptr = matrix; mptr < matrix + DCTSIZE2; mptr++){
    *mptr = *mptr * (*qmatrix);
    qmatrix++;
  }
}

void ChenIDct(int* x, int* y){
  register int i;
  register int* aptr;
  register int a0, a1, a2, a3;
  register int b0, b1, b2, b3;
  register int c0, c1, c2, c3;

  for(i = 0; i < 8; i++){
    aptr = x + i;
    b0 = LS(*aptr, 2);
    aptr += 8;
    a0 = LS(*aptr, 2);
    aptr += 8;
    b2 = LS(*aptr, 2);
    aptr += 8;
    a1 = LS(*aptr, 2);
    aptr += 8;
    b1 = LS(*aptr, 2);
    aptr += 8;
    a2 = LS(*aptr, 2);
    aptr += 8;
    b3 = LS(*aptr, 2);
    aptr += 8;
    a3 = LS(*aptr, 2);

    c0 = MSCALE((c7d16 * a0) - (c1d16 * a3));
    c1 = MSCALE((c3d16 * a2) - (c5d16 * a1));
    c2 = MSCALE((c3d16 * a1) - (c5d16 * a2));
    c3 = MSCALE((c1d16 * a0) - (c7d16 * a3));

    a0 = MSCALE(c1d4 * (b0 + b1));
    a1 = MSCALE(c1d4 * (b0 - b1));
    a2 = MSCALE((c3d8 * b2) - (c1d8 * b3));
    a3 = MSCALE((c1d8 * b2) + (c3d8 * b3));

    b0 = a0 + a3;
    b1 = a1 + a2;
    b2 = a1 - a2;
    b3 = a0 - a3;

    a0 = c0 + c1;
    a1 = c0 - c1;
    a2 = c3 - c2;
    a3 = c3 + c2;

    c0 = a0;
    c1 = MSCALE(c1d4 * (a2 - a1));
    c2 = MSCALE(c1d4 * (a2 + a1));
    c3 = a3;

    aptr = y + i;
    *aptr = b0 + c3;
    aptr += 8;
    *aptr = b1 + c2;
    aptr += 8;
    *aptr = b2 + c1;
    aptr += 8;
    *aptr = b3 + c0;
    aptr += 8;
    *aptr = b3 - c0;
    aptr += 8;
    *aptr = b2 - c1;
    aptr += 8;
    *aptr = b1 - c2;
    aptr += 8;
    *aptr = b0 - c3;
  }

  for(i = 0; i < 8; i++){
    aptr = y + LS(i ,3);
    b0 = *(aptr++);
    a0 = *(aptr++);
    b2 = *(aptr++);
    a1 = *(aptr++);
    b1 = *(aptr++);
    a2 = *(aptr++);
    b3 = *(aptr++);
    a3 = *(aptr);

    c0 = MSCALE((c7d16 * a0) - (c1d16 * a3));
    c1 = MSCALE((c3d16 * a2) - (c5d16 * a1));
    c2 = MSCALE((c3d16 * a1) + (c5d16 * a2));
    c3 = MSCALE((c1d16 * a0) + (c7d16 * a3));

    a0 = MSCALE(c1d4 * (b0 + b1));
    a1 = MSCALE(c1d4 * (b0 - b1));
    a2 = MSCALE((c3d8 * b2) - (c1d8 * b3));
    a3 = MSCALE((c1d8 * b2) + (c3d8 * b3));

    b0 = a0 + a3;
    b1 = a1 + a2;
    b2 = a1 - a2;
    b3 = a0 - a3;

    a0 = c0 + c1;
    a1 = c0 - c1;
    a2 = c3 - c2;
    a3 = c3 + c2;

    c0 = a0;
    c1 = MSCALE(c1d4 * (a2 - a1));
    c2 = MSCALE(c1d4 * (a2 + a1));
    c3 = a3;

    aptr = y + LS(i, 3);
    *(aptr++) = b0 + c3;
    *(aptr++) = b1 + c2;
    *(aptr++) = b2 + c1;
    *(aptr++) = b3 + c0;
    *(aptr++) = b3 - c0;
    *(aptr++) = b2 - c1;
    *(aptr++) = b1 - c2;
    *(aptr) = b0 - c3;
  } 
  for(i = 0, aptr = y; i < 64; i++, aptr++){
    *aptr = (((*aptr < 0) ? (*aptr - 8) : (*aptr + 8)) / 16);
  }
}

void PostshiftIDctMatrix(int* matrix, int shift){
  int* mptr;
  for(mptr = matrix; mptr < matrix + DCTSIZE2; mptr++){
    *mptr += shift;
  }
}

void WriteOneBlock(int* store, unsigned char* out_buf, int width, int height, int voffs, int hoffs){
  int i, e;
  for(i = voffs; i < voffs + DCTSIZE; i++){
    if(i < height){
      int diff;
      diff = width * i;
      for(e = hoffs; e < hoffs + DCTSIZE; e++){
	if(e < width){
	  //printf("%x\n", (*store));
	  out_buf[diff + e] = (unsigned char) (*(store++));
	}
	else{
	  break;
	}
      }
    }
    else{
      break;
    }
  }
}

int DecodeHuffman(int* Xhuff_huffval, int Dhuff_ml, int* Dhuff_maxcode, int* Dhuff_mincode, int* Dhuff_valptr){
  int code, l, p;
  code = buf_getb();
  for(l = 1; code > Dhuff_maxcode[l]; l++){
    code = (code << 1) + buf_getb();
  }
  if(code < Dhuff_maxcode[Dhuff_ml]){
    p = Dhuff_valptr[l] + code - Dhuff_mincode[l];
    return (Xhuff_huffval[p]);
  }
  /*else{
    main_result++;
    printf("Huffman wead error\n");
    exit(4);
    }*/
}

int buf_getv(int n){
  int p, rv;
  n--;
  p = n - read_position;
  while(p > 0){
    if(read_position > 23){
      rv = (current_read_byte << p);
      current_read_byte = pgetc();
      rv |= (current_read_byte >> (8 - p));
      read_position = 7 - p;
      return (rv & lmask[n]);
    }
    current_read_byte = (current_read_byte << 8) | pgetc();
    read_position += 8;
    p -= 8;
  }
  if(!p){
    read_position = -1;
    return (current_read_byte & lmask[n]);
  }
  p = -p;
  read_position = p - 1;
  return ((current_read_byte >> p) & lmask[n]);
}

int buf_getb(){
  if(read_position < 0){
    current_read_byte = pgetc();
    read_position = 7;
  }
  if(current_read_byte & bit_set_mask[read_position--]){
    return (1);
  }
  return (0);
}

static int pgetc(){
  int temp;
  if((temp = *CurHuffReadBuf++) == MARKER_MARKER){
    if((temp = *CurHuffReadBuf++)){
      printf("Unanticipated marker detected.\n");
    }
    else{
      return (MARKER_MARKER);
    }
  }
  return (temp);
}

void BoundIDctMatrix(int* matrix, int Bound){
  int *mptr;
  for(mptr = matrix; mptr < matrix + DCTSIZE2; mptr++){
    if(*mptr < 0){
      *mptr = 0;
    }
    else if(*mptr > Bound){
      *mptr = Bound;
    }
  }
}
