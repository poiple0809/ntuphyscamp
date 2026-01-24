import streamlit as st
import numpy as np
from scipy.integrate import quad

def tunnelling_calculator(E_earned, potential_shape, V0_ev, mul_h, mul_m, prop_list):
  # ----------------- 1. 核心計算函式 -----------------
  def wkb_tunneling_probability(E_joules, m_kg, V_potential_func, x_start, x_end, h_joules_sec=6.62607015e-34):
      """
      使用 WKB 近似法計算穿隧機率 (透射係數 T)。

      參數:
      E_joules (float): 粒子的能量 (單位: 焦耳, J)。
      m_kg (float): 粒子的質量 (單位: 公斤, kg)。
      V_potential_func (function): 位能牆函式 V(x)，必須接受一個參數 x (位置)。
      x_start (float): 積分的起始位置 x1 (位能 V(x1) = E)。
      x_end (float): 積分的結束位置 x2 (位能 V(x2) = E)。
      h_joules_sec (float, optional): 普朗克常數 h。預設為 CODATA 值。

      回傳:
      float: 穿隧機率 T (透射係數)。
      """

      # 設置約化普朗克常數 h-bar
      hbar = h_joules_sec / (2 * np.pi)

      # 檢查穿隧條件：必須確保 E < V(x) 在 (x_start, x_end) 區間內
      # 這裡我們信任使用者輸入的 x_start, x_end 是正確的古典轉折點

      # 被積分函式 (kappa(x) 的兩倍)
      # 被積分函式是 2 * kappa(x)
      def integrand(x):
          V_x = V_potential_func(x)
          # 確保 V_x > E，否則積分無效 (將其設為 0 以避免複數，但在實際應用中，x_start/x_end 應該避免此情況)
          if V_x <= E_joules:
              # 在古典允許區 (V < E) 或古典轉折點 (V = E) 處，積分項為 0
              return 0.0
          else:
              # 這是 2 * kappa(x)
              return (2 * np.sqrt(2 * m_kg * (V_x - E_joules))) / hbar

      # 執行數值積分來計算指數項的絕對值 (K_wkb = 2 * integral(kappa(x) dx))
      # quad 回傳 (積分結果, 誤差估計)
      K_wkb, error = quad(integrand, x_start, x_end)

      # 穿隧機率 T = exp(-K_wkb)
      T = np.exp(-K_wkb)

      return T, K_wkb, error


  # ----------------- 2. 定義位能牆函式 (範例) -----------------

  # --- 矩形位壘 (Rectangular Barrier) ---
  def potential_rectangular(x, V0, L):
      """V(x) = V0 for 0 < x < L, else 0"""
      if 0 <= x <= L:
          return V0
      else:
          return 0.0


  # --- 三角形位壘 (Symmetric Triangular Barrier) ---
  # 假設位壘在 0 到 L 之間，峰值 V_peak 在 L/2 處
  def potential_triangular(x, V_peak, L):
      """V(x) = V_peak * (1 - |2x/L - 1|) for 0 < x < L, else 0"""
      if 0 <= x <= L:
          midpoint = L / 2
          # 線性上升 (0 <= x <= L/2)
          if x <= midpoint:
              return V_peak * (2 * x / L)
          # 線性下降 (L/2 < x <= L)
          else:
              return V_peak * (2 * (L - x) / L)
      else:
          return 0.0


  # --- 拋物線位壘 (Parabolic Barrier, 峰值 V_peak_quad 在 L/2 處) ---
  def potential_quadratic(x, V_peak_quad, L):
      """V(x) = V_peak_quad * (1 - 4*(x - L/2)^2 / L^2) for 0 < x < L, else 0"""
      if 0 <= x <= L:
          # 位能函數: V(x) = V_peak_quad * (1 - 4/L^2 * (x - L/2)^2)
          return V_peak_quad * (1 - (4 / L**2) * (x - L/2)**2)
      else:
          return 0.0


  # ----------------- 3. 範例使用與測試 -----------------
  # 計算剩餘能量
  E_consume = 0
  prop_used = True
  if potential_shape - prop_list[0] == 1:
    E_consume += 40
  elif potential_shape - prop_list[0] == 2:
    E_consume += 60
  elif potential_shape - prop_list[0] == 3:
    E_consume += 150
  elif potential_shape - prop_list[0] == 4:
    E_consume += 160
  elif potential_shape - prop_list[0] == 0:
    E_consume += 0
  else:
    prop_used = False

  if V0_ev + prop_list[1] * 100 == 600:
    E_consume += 0
  elif V0_ev + prop_list[1] * 100 == 500:
    E_consume += 60
  elif V0_ev + prop_list[1] * 100 == 400:
    E_consume += 100
  elif V0_ev + prop_list[1] * 100 == 300:
    E_consume += 200
  else:
    prop_used = False

  if mul_h - prop_list[2] == 1:
    E_consume += 0
  elif mul_h - prop_list[2] == 2:
    E_consume += 120
  elif mul_h - prop_list[2] == 3:
    E_consume += 180
  elif mul_h - prop_list[2] == 4:
    E_consume += 270
  else:
    prop_used = False

  if (mul_m == 1 and prop_list[3] == 0) or (mul_m == 0.25 and prop_list[3]):
    E_consume += 0
  elif (mul_m == 0.25) or (mul_m == 0.1 and prop_list[3]):
    E_consume += 120
  elif (mul_m == 0.1) or (mul_m == 0.063 and prop_list[3]):
    E_consume += 180
  elif mul_m == 0.063:
    E_consume += 270
  else:
    prop_used = False

  E_ev = E_earned - E_consume

  if E_ev <= 0:
    st.error("你破產了！")
    return

  # --- 定義物理常數 ---
  # 粒子的質量 (電子質量)
  m_e = 9.109e-31 * mul_m  # kg
  # 普朗克常數 (預設已在函式中，這裡只是列出)
  h = 6.626e-34 * mul_h  # J·s

  # --- 定義位壘參數 ---
  # 矩形位壘高度 (V0)
  V0_joules = V0_ev * 1.602e-19  # J (將 eV 轉為 J)

  # 位壘寬度
  L_nm = 0.02  # nm
  L_m = L_nm * 1e-9  # m

  # --- 定義粒子能量 ---
  # 粒子能量 E (E < V0 for tunneling)
  E_joules = E_ev * 1.602e-19  # J

  if potential_shape == 0:
    # -----------------------------------------------------------
    # --- 矩形位壘 ---
    # -----------------------------------------------------------

    # 對於矩形位壘，古典轉折點是 0 和 L (因為 E < V0)
    x1_rect = 0.0
    x2_rect = L_m
    V_rect = lambda x: potential_rectangular(x, V0_joules, L_m)

    T_rect, K_rect, error_rect = wkb_tunneling_probability(
        E_joules, m_e, V_rect, x1_rect, x2_rect, h
    )

    #print(f"矩形位壘參數：V0 = {V0_ev} eV, L = {L_nm} nm, E = {E_ev} eV")
    #print(f"WKB 指數項 K: {K_rect:.3e}")
    st.success(f"**穿隧機率 T: {T_rect * 100:.2f} %**")

  elif potential_shape == 2:
    # -----------------------------------------------------------
    # --- 三角形位壘 (假設面積相同，V_peak = 2 * V0) ---
    # -----------------------------------------------------------
    #print("--- 三角形位壘 ---")
    V_peak_joules = 2 * V0_joules

    # 計算三角形位壘的古典轉折點 x1 和 x2
    # V(x) = E
    # 左側斜坡: V_peak * (2x/L) = E => x = E * L / (2 * V_peak)
    x1_tri = E_joules * L_m / (2 * V_peak_joules)
    # 右側斜坡: V_peak * (2(L-x)/L) = E => 2(L-x)/L = E/V_peak => L-x = E*L/(2*V_peak) => x = L - (E*L/(2*V_peak))
    x2_tri = L_m - x1_tri

    V_tri = lambda x: potential_triangular(x, V_peak_joules, L_m)

    T_tri, K_tri, error_tri = wkb_tunneling_probability(
        E_joules, m_e, V_tri, x1_tri, x2_tri, h
    )

    #print(f"三角形位壘參數：V_peak = {V_peak_ev} eV, L = {L_nm} nm, E = {E_ev} eV")
    #print(f"古典轉折點範圍: [{x1_tri:.2e} m, {x2_tri:.2e} m]")
    #print(f"WKB 指數項 K: {K_tri:.3e}")
    st.success(f"**穿隧機率 T: {T_tri * 100:.2f} %**")

  elif potential_shape == 1:
    # -----------------------------------------------------------
    # --- 拋物線位壘 (V_peak_quad = 1.5 V0) ---
    # -----------------------------------------------------------
    V_peak_quad_joules = 1.5 * V0_joules
    if 1 - E_joules / V_peak_quad_joules <= 0:
        T_quad = 1
        st.success(f"**穿隧機率 T: {T_quad * 100:.2f} %**")
        if prop_used == False:
          st.error('別忘了自己有道具!')
        return

    # 計算古典轉折點 x1 和 x2 (V(x) = E)
    # E = V_peak_quad * (1 - 4/L^2 * (x - L/2)^2)
    # (x - L/2)^2 = (L^2/4) * (1 - E/V_peak_quad)
    delta_x = (L_m / 2) * np.sqrt(1 - E_joules / V_peak_quad_joules)
    x1_quad = L_m / 2 - delta_x
    x2_quad = L_m / 2 + delta_x

    V_quad = lambda x: potential_quadratic(x, V_peak_quad_joules, L_m)

    T_quad, K_quad, _ = wkb_tunneling_probability(E_joules, m_e, V_quad, x1_quad, x2_quad, h)

    #print(f"古典轉折點範圍: [{x1_quad:.2e} m, {x2_quad:.2e} m]")
    #print(f"WKB 指數項 K: {K_quad:.3e}")
    st.success(f"**穿隧機率 T: {T_quad * 100:.2f} %**")

  else:
    # -----------------------------------------------------------
    # --- 測試 4: delta function ---
    # -----------------------------------------------------------
    h_bar = h / (2 * np.pi)
    T_delta = 1 / (1 + m_e * (V0_joules * 0.02 * 1e-9) ** 2 / 2 / h_bar ** 2 / E_joules)
    st.success(f"**穿隧機率 T: {T_delta * 100:.2f} %**")

  if prop_used == False:
    st.error('別忘了自己有道具!')

  return

# --- Streamlit 介面設計 ---
st.title('2026台大物理營 WKB近似計算機')
st.subheader("賺到的能量 (單位：eV)")
# 1. 輸入 E_earned
E_earned = st.number_input('E_input', value=300, label_visibility="collapsed")

# 問題決定 prop_list (道具擁有情形)
st.subheader("道具獲得情形")
p1 = st.checkbox("御坂位能塑形器 (電池關道具)")
p2 = st.checkbox("位能高度調整器 (彈簧關道具)")
p3 = st.checkbox("光子能量調整器 (光學關道具)")
p4 = st.checkbox("粒子質量調整器 (質量關道具)")
prop_list = [int(p1), int(p2), int(p3), int(p4)]

st.subheader('商店區')

# 問題決定 potential_shape
st.write('1. 請選擇位能牆的形狀：')
shape_choice = st.radio("shape_choosing",
    options=[0, 1, 2, 3, 4],
    format_func=lambda x: ["矩形(rectangle) (預設)", "拋物線 (parabola) (40eV)", "三角形(triangle) (60eV)", "正δ函數 (150eV)", "負δ函數 (160eV)"][x],
    label_visibility="collapsed")

# 問題決定 V0_ev
st.write('2. 請選擇位能牆底下的總面積：\n(單位：eV * 0.02nm)')
v0_display_labels = {
    600: "600 (預設)",
    500: "500 (60eV)",
    400: "400 (100eV)",
    300: "300 (200eV)"
}

v0_choice = st.radio(
    "v0_choosing",
    options=[600, 500, 400, 300],
    format_func=lambda x: v0_display_labels.get(x),
    label_visibility="collapsed"
)

# 問題決定 mul_h
st.write('3. 請選擇普朗克常數倍數：')
h_display_labels = {
    1: "1 (預設)",
    2: "2 (120eV)",
    3: "3 (180eV)",
    4: "4 (270eV)"
}

h_choice = st.radio('h_choosing',
                    options=[1, 2, 3, 4],
                    format_func=lambda x: h_display_labels.get(x),
                    label_visibility="collapsed")

# 問題決定 mul_m
st.write('4. 請選擇電子質量倍數：')
m_display_labels = {
    1: "1     (預設)",
    0.25: "0.25  (120eV)",
    0.1: "0.1   (180eV)",
    0.063: "0.063 (270eV)"
}

m_choice = st.radio("m_choosing", options=[1, 0.25, 0.1, 0.063],
                    format_func=lambda x: m_display_labels.get(x),
                    label_visibility="collapsed")

# 執行計算
if st.button("計算"):
    # 調用你的函式 (記得修改函式內的 print 為 st.write)
    tunnelling_calculator(E_earned, shape_choice, v0_choice, h_choice, m_choice, prop_list)