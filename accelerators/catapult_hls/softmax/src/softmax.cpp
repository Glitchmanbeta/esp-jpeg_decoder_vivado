// Copyright (c) 2011-2020 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#include "softmax.hpp"

#include <ac_math/ac_softmax_pwl.h>

//
// Compute functions
//

template <class T1, unsigned S1, class T2, unsigned S2>
void compute(plm_t<T1,S1> *input, plm_t<T2,S2> *output) {
    ac_math::ac_softmax_pwl(input->data, output->data);
}

//
// Processes
//
#pragma design modulario<sync>
void softmax::config_accelerator() {
    // HLS_DEFINE_PROTOCOL("config");
    done.write(false); wait();
    //ESP_REPORT_INFO("start configuration");
    // Wait for the configuration signal
    bool end = false;

#pragma hls_unroll no
CONFIG_LOOP:
    do
    {
        wait();
        end = conf_done.read();
    } while (!end);

    // Configuration completed
    done.write(true);

    //ESP_REPORT_INFO("end configuration");

#pragma hls_unroll no
CONFIG_DONE_LOOP:
    while (true) { wait(); }
}

void softmax::load_input() {

    // Load-process reset
    {
        reset_load_input();
        debug = 0;
        wait();
    }

    // Load-process config
    uint32_t batch;
    {
        wait_for_config(); // config process
        conf_info_t config = conf_info.read();

        batch = config.batch;

        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Load config(): batch = %u", ESP_TO_UINT32(batch));
    }

    uint32_t offset = 0;

    // TODO Disable explicit ping-pong buffering. Does Catapult HLS infer
    // ping-pong buffering on its own?
    //bool ping = true;

    ESP_REPORT_TIME(VON, sc_time_stamp(), "load_input(): LOAD_BATCH_LOOP: batch = %u", ESP_TO_UINT32(batch));
    ESP_REPORT_TIME(VON, sc_time_stamp(), "load_input():    LOAD_DATA_INNER_LOOP = %u", PLM_SIZE);

    // Load-process body
LOAD_BATCH_LOOP:
    for (uint32_t b = 0; b < batch; b++) {

        dma_info_t dma_info(offset, PLM_SIZE, 32);

        offset += PLM_SIZE;

        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Load load(): dma_info.index = %u, dma_info.length = %u, dma_info.size = %llu", ESP_TO_UINT32(dma_info.index), ESP_TO_UINT32(dma_info.length), dma_info.size.to_uint64());

        DMA_WRITE(dma_info, dma_read_ctrl);

        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Load load(): dma_read_ctrl done!");

        plm_t<FPDATA_IN, PLM_SIZE> plm_local;

LOAD_DATA_INNER_LOOP:
        for (uint16_t i = 0; i < PLM_SIZE; i++) {

            FPDATA_IN data;
            sc_dt::sc_bv<64> data_bv;
            ac_int<32> data_ac;

            DMA_READ(data_bv, dma_read_chnl);

            // DMA_WIDTH = 64
            // discard bits in the range(63,32)
            // keep bits in the range(31,0)
            data_ac = ac_int<32>(data_bv.range(31,0).to_uint());
            data.set_slc(0, data_ac);
            plm_local.data[i] = data;
        }


        //if (ping) {
        //    plm0_in.write(plm_local);
        //} else {
        //    plm1_in.write(plm_local);
        //}
        plm_in.write(plm_local);
        
        load_compute_handshake();
        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Load load() --> compute()");
        //ping = !ping;
    }

    // Load-process done
    {
        process_done();
    }
}

void softmax::compute_kernel() {

    // Compute-process reset
    {
        reset_compute_kernel();
        wait();
    }

    // Compute-process config
    uint32_t batch;
    {
        wait_for_config(); // config process
        conf_info_t config = conf_info.read();

        batch = config.batch;

        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Compute config(): batch = %u", ESP_TO_UINT32(batch));
    }

    // TODO Disable explicit ping-pong buffering. Does Catapult HLS infer
    // ping-pong buffering on its own?
    //bool ping = true;

    ESP_REPORT_TIME(VON, sc_time_stamp(), "compute_kernel(): COMPUTE_BATCH_LOOP: batch = %u", ESP_TO_UINT32(batch));

    // Compute-process body
COMPUTE_BATCH_LOOP:
    for (uint32_t b = 0; b < batch; b++) {

        compute_load_handshake();
        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Compute compute() <---> load()");

        plm_t<FPDATA_IN, PLM_SIZE> plm_local_in;
        plm_t<FPDATA_OUT, PLM_SIZE> plm_local_out;

        //if (ping) {
        //    plm_local_in = plm0_in.read();
        //} else {
        //    plm_local_in = plm1_in.read();
        //}
        plm_local_in = plm_in.read();

        compute<FPDATA_IN, PLM_SIZE, FPDATA_OUT, PLM_SIZE>(&plm_local_in, &plm_local_out);

        //if (ping) {
        //    plm0_out.write(plm_local_out);
        //} else {
        //    plm1_out.write(plm_local_out);
        //}
        plm_out.write(plm_local_out);

        compute_store_handshake();
        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Compute compute() ---> store()");

        //ping = !ping;
    }

    // Compute-process done
    {
        process_done();
    }
}

void softmax::store_output() {

    // Store-process reset
    {
        reset_store_output();
        wait();
    }

    // Store-process config
    uint32_t batch;
    {
        wait_for_config(); // config process
        conf_info_t config = conf_info.read();

        batch = config.batch;

        ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Store config(): batch = %u", ESP_TO_UINT32(batch));
    }

    uint32_t offset = PLM_SIZE * batch;

    // TODO Disable explicit ping-pong buffering. Does Catapult HLS infer
    // ping-pong buffering on its own?
    //bool ping = true;

    ESP_REPORT_TIME(VON, sc_time_stamp(), "store_output(): STORE_BATCH_LOOP: batch = %u", ESP_TO_UINT32(batch));
    ESP_REPORT_TIME(VON, sc_time_stamp(), "store_output():    STORE_DATA_INNER_LOOP = %u", PLM_SIZE);

    // Store-process body
STORE_BATCH_LOOP:
    for (uint32_t b = 0; b < batch; b++) {

            store_compute_handshake();
            ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Store store() --> compute()");

            dma_info_t dma_info(offset, PLM_SIZE, 32);

            offset += PLM_SIZE;

            ESP_REPORT_TIME(VOFF, sc_time_stamp(), "Store store(): dma_info.index = %u, dma_info.length = %u, dma_info.size = %llu", ESP_TO_UINT32(dma_info.index), ESP_TO_UINT32(dma_info.length), dma_info.size.to_uint64());

            DMA_WRITE(dma_info, dma_write_ctrl);

            plm_t<FPDATA_OUT, PLM_SIZE> plm_local;

            //if (ping) {
            //    plm_local = plm0_out.read();
            //} else {
            //    plm_local = plm1_out.read();
            //}
            plm_local = plm_out.read();

STORE_DATA_INNER_LOOP:
            for (uint16_t i = 0; i < PLM_SIZE; i++) {

                FPDATA_OUT data = plm_local.data[i];

                // DMA_WIDTH = 64
                // set to a constante value range(63,32)
                // return results on the range(31,0)
                sc_dt::sc_bv<64> data_bv;
                data_bv.range(63,32) = sc_dt::sc_bv<32>(0xdeadbeef);
                data_bv.range(31,0) = data.template slc<32>(0);

                DMA_WRITE(data_bv, dma_write_chnl);
            }


            //ping = !ping;
    }

    // Store-process done
    {
        accelerator_done();
        process_done();
    }
}

// ***************************************************
// *** YOU SHOULD NOT EDIT THE FOLLOWING FUNCTIONS ***
// ***************************************************

//
// Reset functions
//

inline void softmax::reset_dma_read() {
    DMA_WRITE_RESET(dma_read_ctrl);
    DMA_READ_RESET(dma_read_chnl);
}

inline void softmax::reset_dma_write() {
    DMA_WRITE_RESET(dma_write_ctrl);
    DMA_WRITE_RESET(dma_write_chnl);
}

inline void softmax::reset_accelerator_done() {
    acc_done.write(false);
}

//
// Functions
//

inline void softmax::reset_load_input() {
    input_ready.reset_req();
    reset_dma_read();
}

inline void softmax::reset_compute_kernel() {
    input_ready.reset_ack();
    output_ready.reset_req();
}

inline void softmax::reset_store_output() {
    reset_accelerator_done();
    reset_dma_write();
    output_ready.reset_ack();
}

inline void softmax::load_compute_handshake() {
    input_ready.req();
}

inline void softmax::compute_load_handshake() {
    input_ready.ack();
}

inline void softmax::compute_store_handshake() {
    output_ready.req();
}

inline void softmax::store_compute_handshake() {
    output_ready.ack();
}

inline void softmax::wait_for_config() {
#pragma hls_unroll no
WAIT_FOR_CONFIG_LOOP:
    while (!done.read()) { wait(); }
}

inline void softmax::process_done() {
#pragma hls_unroll no
PROCESS_DONE_LOOP:
    do { wait(); } while (true);
}

inline void softmax::accelerator_done() {
    acc_done.write(true); wait(); wait(); wait();
    acc_done.write(false);
}