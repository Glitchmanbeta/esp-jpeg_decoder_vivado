// Copyright (c) 2011-2021 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include <linux/of_device.h>
#include <linux/mm.h>

#include <asm/io.h>

#include <esp_accelerator.h>
#include <esp.h>

#include "jpeg_decoder_vivado.h"

#define DRV_NAME	"jpeg_decoder_vivado"

/* <<--regs-->> */
#define JPEG_DECODER_TEST_REG 0x40

struct jpeg_decoder_vivado_device {
	struct esp_device esp;
};

static struct esp_driver jpeg_decoder_driver;

static struct of_device_id jpeg_decoder_device_ids[] = {
	{
		.name = "SLD_JPEG_DECODER_VIVADO",
	},
	{
		.name = "eb_056",
	},
	{
		.compatible = "sld,jpeg_decoder_vivado",
	},
	{ },
};

static int jpeg_decoder_devs;

static inline struct jpeg_decoder_vivado_device *to_jpeg_decoder(struct esp_device *esp)
{
	return container_of(esp, struct jpeg_decoder_vivado_device, esp);
}

static void jpeg_decoder_prep_xfer(struct esp_device *esp, void *arg)
{
	struct jpeg_decoder_vivado_access *a = arg;

	/* <<--regs-config-->> */
	iowrite32be(a->test, esp->iomem + JPEG_DECODER_TEST_REG);
	iowrite32be(a->src_offset, esp->iomem + SRC_OFFSET_REG);
	iowrite32be(a->dst_offset, esp->iomem + DST_OFFSET_REG);

}

static bool jpeg_decoder_xfer_input_ok(struct esp_device *esp, void *arg)
{
	/* struct jpeg_decoder_vivado_device *jpeg_decoder = to_jpeg_decoder(esp); */
	/* struct jpeg_decoder_vivado_access *a = arg; */

	return true;
}

static int jpeg_decoder_probe(struct platform_device *pdev)
{
	struct jpeg_decoder_vivado_device *jpeg_decoder;
	struct esp_device *esp;
	int rc;

	jpeg_decoder = kzalloc(sizeof(*jpeg_decoder), GFP_KERNEL);
	if (jpeg_decoder == NULL)
		return -ENOMEM;
	esp = &jpeg_decoder->esp;
	esp->module = THIS_MODULE;
	esp->number = jpeg_decoder_devs;
	esp->driver = &jpeg_decoder_driver;
	rc = esp_device_register(esp, pdev);
	if (rc)
		goto err;

	jpeg_decoder_devs++;
	return 0;
 err:
	kfree(jpeg_decoder);
	return rc;
}

static int __exit jpeg_decoder_remove(struct platform_device *pdev)
{
	struct esp_device *esp = platform_get_drvdata(pdev);
	struct jpeg_decoder_vivado_device *jpeg_decoder = to_jpeg_decoder(esp);

	esp_device_unregister(esp);
	kfree(jpeg_decoder);
	return 0;
}

static struct esp_driver jpeg_decoder_driver = {
	.plat = {
		.probe		= jpeg_decoder_probe,
		.remove		= jpeg_decoder_remove,
		.driver		= {
			.name = DRV_NAME,
			.owner = THIS_MODULE,
			.of_match_table = jpeg_decoder_device_ids,
		},
	},
	.xfer_input_ok	= jpeg_decoder_xfer_input_ok,
	.prep_xfer	= jpeg_decoder_prep_xfer,
	.ioctl_cm	= JPEG_DECODER_VIVADO_IOC_ACCESS,
	.arg_size	= sizeof(struct jpeg_decoder_vivado_access),
};

static int __init jpeg_decoder_init(void)
{
	return esp_driver_register(&jpeg_decoder_driver);
}

static void __exit jpeg_decoder_exit(void)
{
	esp_driver_unregister(&jpeg_decoder_driver);
}

module_init(jpeg_decoder_init)
module_exit(jpeg_decoder_exit)

MODULE_DEVICE_TABLE(of, jpeg_decoder_device_ids);

MODULE_AUTHOR("Emilio G. Cota <cota@braap.org>");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("jpeg_decoder_vivado driver");
