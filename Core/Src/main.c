/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2022 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include <math.h>

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */


#define M 	1	 			// Cart Mass[kg]
#define m 	0.14 			// Pendulum Mass[kg]
#define l 	60	 			// Pendulum length[cm]
#define b   0.5  			// Cart friction constant[kg/s]
#define r   1.27 			// Movement distance per rotation of the ball screw[cm]
#define km  4.9  			// Motor torque constant[Ncm/A]
#define kb  0.0507  		// Motor Back EMF constant[V/rad/s]
#define R   0.3     		// Motor Resistance[ohm]
#define F_v 80.8074			// (2*pi*Km)/(r*R)
#define F_r 20.7691			// b+(2*pi/r)^2*Km*Kb/R
#define I 	168				// Pendulum Inertia

#define pi 	3.141592653589793238462643383279502884197169399375105820974944
#define g 	981
#define dt_vel 0.001		// 1ms 	Position & velocity & accel / Encoder
#define dt_ome 0.02			// 20ms
#define RAD2DEG 57.295779
#define DEG2RAD 0.0174532

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
TIM_HandleTypeDef htim1;
TIM_HandleTypeDef htim2;
TIM_HandleTypeDef htim3;
TIM_HandleTypeDef htim9;

UART_HandleTypeDef huart2;

/* USER CODE BEGIN PV */

// Setting
volatile uint32_t cnt_1ms;
volatile uint32_t cnt_20ms;
volatile int mode = 0;
volatile int flag_LQR = 0;
volatile int flag_Swing = 0;

//////////encoder variable//////////

// Motor Encoder
volatile uint16_t encoder_cnt_m = 0;
volatile uint16_t encoder_cnt_pre_m = 0;
volatile int32_t encoder_err_m = 0;
volatile double angle_m = 0;
volatile double angle_pre_m = 0;
volatile double omega_m = 0;

// Pendulum Encoder
volatile uint16_t encoder_cnt_b = 0;
volatile uint16_t encoder_cnt_pre_b = 0;
volatile int32_t encoder_err_b = 0;
volatile double angle_b = -180;
volatile double angle_pre_b = 0;
volatile double omega_b = 0;
volatile double Motor_counter = 0;

// Cart
volatile double distance_c = 0;
volatile double distance_pre_c = 0;
volatile double velocity_c = 0;
volatile double velocity_pre_c = 0;

// PWM
volatile uint16_t Motor_ccr = 0;

//LQR
// Cart Position   Cart Velocity   Pendulum Position   Pendulum Velocity
//volatile double K[4] = {-15.2361,  -65.8163,  592.5298,  10.9238}; // best
volatile double K[4] = {-1.0000,   -1.1004,  332.2568,  88.7627}; // Hyeon-min
//volatile double K[4] = {-1.0000,   -1.1004,  592.2568,  40.7627}; // best
volatile double LQR = 0;
volatile double angle_gap = 0;
volatile int cnt_LQR = 0;

//Swing up
volatile double Swing = 0;
volatile double G = 0;
volatile double V_pre = 0;
volatile double E_p = 0;
volatile int sign = 0;
volatile double u = 0;
volatile double w = 2.31;
volatile double u_a = 8.5;

//cos
volatile double cos_angle_b = 0;
volatile double max_cos = -2;

//LPF
volatile double LPF_fc_u = 1;
volatile double LPF_fc = 8;
volatile double LPF = 0;
volatile double LPF_tau = 0;
volatile double LPF_past = 0;
volatile double dt = 0.02;

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART2_UART_Init(void);
static void MX_TIM1_Init(void);
static void MX_TIM2_Init(void);
static void MX_TIM3_Init(void);
static void MX_TIM9_Init(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */
void Get_max_cos(void) 	// for Control Switching, Max of cos value
{
	if(cos_angle_b > max_cos) max_cos = cos_angle_b;
}
void Get_Ep(void) 		// for Energy of mechanical value
{
	E_p = ((I * pow(omega_b,2)) / 2) + m * g * l * (cos(angle_b * DEG2RAD) - 1);
}
int sgn(double x)		// determine sign
{
	if(x < 0) sign = -1;
	if(x > 0) sign = 1;

	return sign;
}
double Low_Pass_Filter(double data)
{
	LPF_tau = 1 / (2 * LPF_fc * pi);
	LPF = (dt * data + LPF_tau * LPF_past) / (LPF_tau + dt);
	LPF_past = LPF;

	return LPF;
}
double Low_Pass_Filter_u(double data)
{
	LPF_tau = 1 / (2 * LPF_fc_u * pi);
	LPF = (dt * data + LPF_tau * LPF_past) / (LPF_tau + dt);
	LPF_past = LPF;

	return LPF;
}
void Get_u(void)		// acceleration of the cart = u
{
	u = (((u_a * fabs( (sgn(E_p) * omega_b * cos(angle_b * DEG2RAD)) - 2 * w * velocity_c) ) + (2 * w * distance_c * velocity_c) ) / ((sgn(E_p) * omega_b * cos(angle_b * DEG2RAD)) - 2 * w * velocity_c)) / 6;
    if(((sgn(E_p) * omega_b * cos(angle_b * DEG2RAD)) - (2 * w * velocity_c)) == 0)
	{
    	u = 1;
	}
}
void Swing_up(void)		// for Swing-up Algorithm with Lyapunov function and decrease of Energy of Pendulum
{
	Swing = G * ( (M + (m * pow(sin(angle_b * DEG2RAD), 2))) * u + ( m * g * sin(angle_b * DEG2RAD) * cos(angle_b * DEG2RAD)) + (F_r * velocity_c) - (m * l * sin(angle_b * DEG2RAD) * pow((omega_b), 2)) ) / F_v;
}
void LQR_Control(void)	// for LQR control
{
	if(angle_b > -180)
	{
		LQR = (K[0] * (-distance_c)) + (K[1] * (-velocity_c)) + (K[2] * (-angle_b)) + (K[3] * (omega_b));
	}
	else
	{
		LQR = (K[0] * (-distance_c)) + (K[1] * (-velocity_c)) + (K[2] * (-angle_b - 360)) + (K[3] * (omega_b));
	}
}
void LQR_Control_Stabilize(void)	// for LQR control
{
	//K[4] = {-15.2361,  -65.8163,  592.5298,  10.9238};
	K[0] = -15.2361;
	K[1] = -65.8163;
	K[2] = 632.5298;
	K[3] = 15.9238;
	if(angle_b > -180)
		{
			LQR = (K[0] * (-distance_c)) + (K[1] * (-velocity_c)) + (K[2] * (-angle_b)) + (K[3] * (-omega_b));
		}
		else
		{
			LQR = (K[0] * (-distance_c)) + (K[1] * (-velocity_c)) + (K[2] * (-angle_b - 360)) + (K[3] * (-omega_b));
		}
}
void Make_PWM_LQR(double input)		//
{
	if(input > 0) 	//?��방향
	{
		Motor_ccr = (double)(input);
		HAL_GPIO_WritePin(GPIOC, GPIO_PIN_4, GPIO_PIN_RESET);
		if(Motor_ccr > 8999) Motor_ccr = 8999;
		TIM3->CCR1 = Motor_ccr;
	}
	else 			//?��방향
	{
		Motor_ccr = (double)(-input);
		HAL_GPIO_WritePin(GPIOC, GPIO_PIN_4, GPIO_PIN_SET);
		if(Motor_ccr > 8999) Motor_ccr = 8999;
		TIM3->CCR1 = Motor_ccr;
	}
}
void Make_PWM_Swing(double input)	//
{
	if(input > 0) 	//?��방향
	{
		Motor_ccr = (double)(input);
		HAL_GPIO_WritePin(GPIOC, GPIO_PIN_4, GPIO_PIN_SET);
		if(Motor_ccr > 8999) Motor_ccr = 8999;
		TIM3->CCR1 = Motor_ccr;
	}
	else 			//?��방향
	{
		Motor_ccr = (double)(-input);
		HAL_GPIO_WritePin(GPIOC, GPIO_PIN_4, GPIO_PIN_RESET);
		if(Motor_ccr > 8999) Motor_ccr = 8999;
		TIM3->CCR1 = Motor_ccr;
	}
}
void HAL_GPIO_EXTI_Callback(uint16_t GPIO_Pin) 	// Button interrupt
{
	if(GPIO_Pin == GPIO_PIN_13)
	{
		HAL_GPIO_TogglePin(LD2_GPIO_Port,LD2_Pin);
		if(mode == 0) mode = 1;

		// test
		//HAL_GPIO_TogglePin(GPIOC, GPIO_PIN_4);
		//TIM3->CCR1 = 1000;

	}
}
/////////////////////Control Start//////////////////////////////
////////////////////////////////////////////////////////////////
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim)
{
/////////////////////Encoder Sensing////////////////////////////
////////////////////////////////////////////////////////////////
	if(htim == &htim9 && mode == 1)
	{
		cnt_1ms++;

		// encoder & angle & omega sensing
			// Pendulum
		angle_pre_b = angle_b;
		encoder_cnt_pre_b = encoder_cnt_b;
		encoder_cnt_b = TIM1->CNT;
		encoder_err_b = encoder_cnt_b - encoder_cnt_pre_b;
			// Motor
		angle_pre_m = angle_m;
		encoder_cnt_pre_m = encoder_cnt_m;
		encoder_cnt_m = TIM2->CNT;
		encoder_err_m = encoder_cnt_m - encoder_cnt_pre_m;

			// Cart
		distance_pre_c = distance_c;
		velocity_pre_c = velocity_c;

		// encoder exception handling
			// Pendulum
		if(encoder_cnt_b>65000 && encoder_cnt_pre_b < 10000)  encoder_err_b -= 65535 ;
		else if(encoder_cnt_b< 10000 && encoder_cnt_pre_b > 65000)  encoder_err_b += 65535 ;
		angle_b += (encoder_err_b/ 8000.) * 360; 	// degree
        omega_b = (angle_b - angle_pre_b) / dt_vel; // degree/s

        	// Motor
        if(encoder_cnt_m>60000 && encoder_cnt_pre_m < 10000)  encoder_err_m -= 65535 ;
        else if(encoder_cnt_m < 10000 && encoder_cnt_pre_m > 60000)  encoder_err_m += 65535 ;
        angle_m += (encoder_err_m/ 8000.) * 360;	// degree
        omega_m = (angle_m - angle_pre_m) / dt_vel; // degre/s*/

        	// Cart
        Motor_counter = ((angle_m / 360) ); 				 	// Number of Revolutions
        distance_c = Motor_counter * 1.27;			 			// Cart distance cm
        velocity_c = (distance_c - distance_pre_c) / dt_vel; 	// Cart velocity cm/s

        	// Swing
        cos_angle_b = cos(angle_b * DEG2RAD);	//cos
        Get_max_cos();

///////////////////////Control Period/////////////////////////////
//////////////////////////////////////////////////////////////////
		if (cnt_1ms % 20 == 0)
		{
			cnt_20ms++;

//////////////////////////////////////////////////////////////////
			// Swing-up
			if (max_cos < 0.96)flag_Swing = 1;			// +- 25
			if (flag_Swing == 1 && cos_angle_b < 0.94)	// +- 25
			{
				if (cnt_20ms < 20)
				{
					Make_PWM_Swing(-6500);
				}
				else if (cnt_20ms > 20 && cnt_20ms < 45)
				{
					Make_PWM_Swing(6500);
				}
				else
				{
					Get_Ep();
					Get_u();
					u = Low_Pass_Filter_u(u);
					// 1
					if (max_cos < 0.) G = 0.5;
					// 2
					if (max_cos > 0. && max_cos < 0.6) G = 0.2;
					// 3
					if (max_cos > 0.6) G = 0.063;
					// Cala input
					Swing_up();
					Swing = Low_Pass_Filter(Swing);
					Make_PWM_Swing(Swing);
				}
			}
//////////////////////////////////////////////////////////////////
			// LQR
			if (max_cos > 0.96) flag_LQR = 1;			// +- 25
			if (flag_LQR == 1 && cos_angle_b > 0.94) 	// +- 25
			{
				cnt_LQR++;
				flag_Swing = 0;
				Swing = 0;
				if(cnt_LQR < 8)
				{
					LQR_Control();
					Make_PWM_LQR(LQR);
				}
				else
				{
					LQR_Control_Stabilize();
					Make_PWM_LQR(LQR);
				}
			}
//////////////////////////////////////////////////////////////////
			// LQR --- Swing-up
			if(max_cos > 0.98 && cos_angle_b < 0.5)
			{
				flag_LQR = 0;
				LQR = 0;
				flag_Swing = 1;
			 	Make_PWM_Swing(0);
			}
		}
	}
}
/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_USART2_UART_Init();
  MX_TIM1_Init();
  MX_TIM2_Init();
  MX_TIM3_Init();
  MX_TIM9_Init();
  /* USER CODE BEGIN 2 */
  HAL_TIM_Base_Start_IT(&htim9);
  HAL_TIM_Encoder_Start(&htim1, TIM_CHANNEL_ALL);	//bong 	= Timer1
  HAL_TIM_Encoder_Start(&htim2, TIM_CHANNEL_ALL);	//motor = Timer2
  HAL_TIM_PWM_Start(&htim3,TIM_CHANNEL_1);
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
   while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 8;
  RCC_OscInitStruct.PLL.PLLN = 180;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 2;
  RCC_OscInitStruct.PLL.PLLR = 2;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Activate the Over-Drive mode
  */
  if (HAL_PWREx_EnableOverDrive() != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV4;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV2;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_5) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief TIM1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM1_Init(void)
{

  /* USER CODE BEGIN TIM1_Init 0 */

  /* USER CODE END TIM1_Init 0 */

  TIM_Encoder_InitTypeDef sConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM1_Init 1 */

  /* USER CODE END TIM1_Init 1 */
  htim1.Instance = TIM1;
  htim1.Init.Prescaler = 0;
  htim1.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim1.Init.Period = 65535;
  htim1.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim1.Init.RepetitionCounter = 0;
  htim1.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  sConfig.EncoderMode = TIM_ENCODERMODE_TI12;
  sConfig.IC1Polarity = TIM_ICPOLARITY_RISING;
  sConfig.IC1Selection = TIM_ICSELECTION_DIRECTTI;
  sConfig.IC1Prescaler = TIM_ICPSC_DIV1;
  sConfig.IC1Filter = 0;
  sConfig.IC2Polarity = TIM_ICPOLARITY_RISING;
  sConfig.IC2Selection = TIM_ICSELECTION_DIRECTTI;
  sConfig.IC2Prescaler = TIM_ICPSC_DIV1;
  sConfig.IC2Filter = 0;
  if (HAL_TIM_Encoder_Init(&htim1, &sConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim1, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM1_Init 2 */

  /* USER CODE END TIM1_Init 2 */

}

/**
  * @brief TIM2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM2_Init(void)
{

  /* USER CODE BEGIN TIM2_Init 0 */

  /* USER CODE END TIM2_Init 0 */

  TIM_Encoder_InitTypeDef sConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};

  /* USER CODE BEGIN TIM2_Init 1 */

  /* USER CODE END TIM2_Init 1 */
  htim2.Instance = TIM2;
  htim2.Init.Prescaler = 0;
  htim2.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim2.Init.Period = 65535;
  htim2.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim2.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  sConfig.EncoderMode = TIM_ENCODERMODE_TI12;
  sConfig.IC1Polarity = TIM_ICPOLARITY_RISING;
  sConfig.IC1Selection = TIM_ICSELECTION_DIRECTTI;
  sConfig.IC1Prescaler = TIM_ICPSC_DIV1;
  sConfig.IC1Filter = 0;
  sConfig.IC2Polarity = TIM_ICPOLARITY_RISING;
  sConfig.IC2Selection = TIM_ICSELECTION_DIRECTTI;
  sConfig.IC2Prescaler = TIM_ICPSC_DIV1;
  sConfig.IC2Filter = 0;
  if (HAL_TIM_Encoder_Init(&htim2, &sConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim2, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM2_Init 2 */

  /* USER CODE END TIM2_Init 2 */

}

/**
  * @brief TIM3 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM3_Init(void)
{

  /* USER CODE BEGIN TIM3_Init 0 */

  /* USER CODE END TIM3_Init 0 */

  TIM_MasterConfigTypeDef sMasterConfig = {0};
  TIM_OC_InitTypeDef sConfigOC = {0};

  /* USER CODE BEGIN TIM3_Init 1 */

  /* USER CODE END TIM3_Init 1 */
  htim3.Instance = TIM3;
  htim3.Init.Prescaler = 1;
  htim3.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim3.Init.Period = 4199;
  htim3.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim3.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_PWM_Init(&htim3) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_RESET;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim3, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sConfigOC.OCMode = TIM_OCMODE_PWM1;
  sConfigOC.Pulse = 0;
  sConfigOC.OCPolarity = TIM_OCPOLARITY_HIGH;
  sConfigOC.OCFastMode = TIM_OCFAST_DISABLE;
  if (HAL_TIM_PWM_ConfigChannel(&htim3, &sConfigOC, TIM_CHANNEL_1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM3_Init 2 */

  /* USER CODE END TIM3_Init 2 */
  HAL_TIM_MspPostInit(&htim3);

}

/**
  * @brief TIM9 Initialization Function
  * @param None
  * @retval None
  */
static void MX_TIM9_Init(void)
{

  /* USER CODE BEGIN TIM9_Init 0 */

  /* USER CODE END TIM9_Init 0 */

  TIM_ClockConfigTypeDef sClockSourceConfig = {0};

  /* USER CODE BEGIN TIM9_Init 1 */

  /* USER CODE END TIM9_Init 1 */
  htim9.Instance = TIM9;
  htim9.Init.Prescaler = 180-1;
  htim9.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim9.Init.Period = 1000-1;
  htim9.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim9.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim9) != HAL_OK)
  {
    Error_Handler();
  }
  sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
  if (HAL_TIM_ConfigClockSource(&htim9, &sClockSourceConfig) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN TIM9_Init 2 */

  /* USER CODE END TIM9_Init 2 */

}

/**
  * @brief USART2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART2_UART_Init(void)
{

  /* USER CODE BEGIN USART2_Init 0 */

  /* USER CODE END USART2_Init 0 */

  /* USER CODE BEGIN USART2_Init 1 */

  /* USER CODE END USART2_Init 1 */
  huart2.Instance = USART2;
  huart2.Init.BaudRate = 115200;
  huart2.Init.WordLength = UART_WORDLENGTH_8B;
  huart2.Init.StopBits = UART_STOPBITS_1;
  huart2.Init.Parity = UART_PARITY_NONE;
  huart2.Init.Mode = UART_MODE_TX_RX;
  huart2.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart2.Init.OverSampling = UART_OVERSAMPLING_16;
  if (HAL_UART_Init(&huart2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART2_Init 2 */

  /* USER CODE END USART2_Init 2 */

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOC, GPIO_PIN_4, GPIO_PIN_RESET);

  /*Configure GPIO pin : PC13 */
  GPIO_InitStruct.Pin = GPIO_PIN_13;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /*Configure GPIO pin : LD2_Pin */
  GPIO_InitStruct.Pin = LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LD2_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : PC4 */
  GPIO_InitStruct.Pin = GPIO_PIN_4;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /* EXTI interrupt init*/
  HAL_NVIC_SetPriority(EXTI15_10_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);

}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
