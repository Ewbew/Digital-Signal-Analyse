



/* tmc2130_driver.c */
struct spi_device *my_spi_device; // Assume already initialized

/* Attribute(s) */
/* YOUR CODE HERE */

static ssize_t over_temperature_show(struct device *dev, struct device_attribute *attr, char *buf){
    u8 o_temp_state = 0;
    int err;

    tmc2130_spi_read_ot(dev, &o_temp_state);

    if(o_temp_state){
        int len = sprintf(buf, "%d\n", 1);
    } else {
        int len = sprintf(buf, "%d\n", 0);
    }
    return err;
}

DEVICE_ATTR_RO(over_temperature);

static struct attribute *tmc2130_attrs[] = {
   &dev_attr_over_temperature.attr,
   NULL,
};

ATTRIBUTE_GROUPS(tmc_2130);

////
    /* Retrieve drvdata ptr (set in device_create) */
    struct gpio_dev *d = dev_get_drvdata(dev);
    int value = d->toggle_state;
    int len = sprintf(buf, "%d\n", value);
  
    return len;





static struct device_attribute dev_attr_led_state = {
.attr = {.name = "led_state",
.mode = 0644},
.show = led_state_show,
.store = led_state_store,
};



static ssize_t my_gpio_state_show(struct device *dev, struct device_attribute *attr, char *buf) 
{ 
    /* Retrieve drvdata ptr (set in device_create) */
    struct gpio_dev *d = dev_get_drvdata(dev);
    int value = d->toggle_state;
    int len = sprintf(buf, "%d\n", value);
  
    return len;
}


























// HAL eksamen s2020
int tmc2130_spi_read_ot(struct spi_device* spi, u8* ot) {
  int err = 0;
  struct spi_transfer t[6]; // I alt 6 bytes tranfer, inklusive startbit og addresse bits
  struct spi_message m;

  	memset(t,0, sizeof(t)); /* Init Memory */
	spi_message_init(&m);	/* Init Msg */
	m.spi = spi;			/* Use current SPI I/F */

    u8 read_bit = 0x01;
    u8 adr = 0xaa; // Address unknown
    u8 data0, data1, data2, data3;

    /* Add first message to queue */
	t[0].tx_buf = &read_bit;
	t[0].rx_buf = NULL;
	t[0].len = 1;
	spi_message_add_tail(&t[0], &m); /* Add message to queue */

    /* Add second message to queue */
	t[1].tx_buf = &adr;
	t[1].rx_buf = NULL;
	t[1].len = 1;
	spi_message_add_tail(&t[1], &m); /* Add message to queue */

    /* Add third message to queue */
	t[2].tx_buf = NULL;
	t[2].rx_buf = &data0;
	t[2].len = 1;
	spi_message_add_tail(&t[2], &m); /* Add message to queue */

    /* Add fourth message to queue */
	t[3].tx_buf = NULL;
	t[3].rx_buf = &data1;
	t[3].len = 1;
	spi_message_add_tail(&t[3], &m); /* Add message to queue */

    /* Add fifth message to queue */
	t[4].tx_buf = NULL;
	t[4].rx_buf = &data2;
	t[4].len = 1;
	spi_message_add_tail(&t[4], &m); /* Add message to queue */

    /* Add sixth message to queue */
	t[5].tx_buf = NULL;
	t[5].rx_buf = &data3;
	t[5].len = 1;
	spi_message_add_tail(&t[5], &m); /* Add message to queue */

    /* Transmit SPI Data */
	err = spi_sync(m.spi, &m);

    /* Only need bit 25, which is in data0, since MSB is transmitted first.
    data0 consists of bits [31:24], so the ot bit must be bit 1 in this u8: */
    *ot = (data0 & 0b00000010);
    /* Bit 1 in ot is therefore now high, if the Over Temperature alarm (the ot bit on the TCM2130) was set during reading*/

  return err;
}


int ReadADC(struct spi_device *spi, int channel, u16 *result){
	struct spi_transfer t[3]; /* Three-byte transfer */
	struct spi_message m;
	int err;

	/* Fejlhåndtering, hvis channel-parameter er forskellig fra 0 el. 1 */
	if (channel != 0 && channel != 1) {
		printk(KERN_ERR "spi_drv: Invalid ADC channel %d\n", channel);
		return -EINVAL;
	}	

	memset(t,0, sizeof(t)); /* Init Memory */
	spi_message_init(&m);	/* Init Msg */
	m.spi = spi;			/* Use current SPI I/F */

	/* Master transmission bytes */
	u8 empty = 0b0000000;	/* Empty data transmission */
	u8 mode = 0b10000000 | (channel << 6);	/* Mode-command transmission */
	u8 start = 0b0000001;	/* Start-command transmission */

	/* Master receiving bytes */
	u8 data1;	/* First 4-bit ADC-value of 12-bit from MCP2032 */
	u8 data2;	/* Last 8-bit ADC-value of 12-bit from MCP2032 */

	/* Add first message to queue */
	t[0].tx_buf = &start;
	t[0].rx_buf = &empty;
	t[0].len = 1;
	spi_message_add_tail(&t[0], &m); /* Add message to queue */

	/* Add second message to queue */
	t[1].tx_buf = &mode;
	t[1].rx_buf = &data1;
	t[1].len = 1;
	spi_message_add_tail(&t[1], &m); /* Add message to queue */

	/* Add third (last) message to queue */
	t[2].tx_buf = &empty;
	t[2].rx_buf = &data2;
	t[2].len = 1;
	spi_message_add_tail(&t[2], &m); /* Add message to queue */

	/* Transmit SPI Data (blocking) */
	err = spi_sync(m.spi, &m);

	/* Decode two data-bytes into 12-bit result and convert to mV */
	*result = ( (data1 & 0b00001111) << 8 ) | data2;
	*result = (*result*MCP3202_VREF_MV)/MCP3202_ADC_MAX;

	return err;}





// First exercise
    /{  compatible = "brcm,bcm2835", "brcm,bcm2836", "brcm,bcm2708", "brcm,bcm2709";  fragment@1 {    target = <&spi0>;    __overlay__ {      #address-cells = <1>;      #size-cells = <0>;      tmc2130: tmc2130@0 {        compatible = "au-ece, tmc2130";        reg = <0>; // spi0.0
        spi-cpha; // CPHA flag is set high (CPHA = 1)        spi-cpol; // CPOL flag is set high (CPHA = 1)        // When cpha and cpol both are high, the SPI mode is 3           spi-max-frequency = <4000000>; // Max frequency using internal clock, 4 MHz      };    };  };}; 