// Copyright (c) 2011-2021 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#ifndef __ESP_CFG_000_H__
#define __ESP_CFG_000_H__

#include "libesp.h"
#include "jpeg_decoder_vivado.h"

typedef int64_t token_t;

/* <<--params-def-->> */
#define TEST 42

/* <<--params-->> */
const int32_t test = TEST;

#define NACC 1

struct jpeg_decoder_vivado_access jpeg_decoder_cfg_000[] = {
	{
		/* <<--descriptor-->> */
		.test = TEST,
		.src_offset = 0,
		.dst_offset = 0,
		.esp.coherence = ACC_COH_NONE,
		.esp.p2p_store = 0,
		.esp.p2p_nsrcs = 0,
		.esp.p2p_srcs = {"", "", "", ""},
	}
};

esp_thread_info_t cfg_000[] = {
	{
		.run = true,
		.devname = "jpeg_decoder_vivado.0",
		.ioctl_req = JPEG_DECODER_VIVADO_IOC_ACCESS,
		.esp_desc = &(jpeg_decoder_cfg_000[0].esp),
	}
};

#endif /* __ESP_CFG_000_H__ */
