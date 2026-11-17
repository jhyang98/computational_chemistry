# soft_east — (s, g) 편향 softened East 모델 계산 도구

활동도 편향 s와 에너지 편향 g를 건 softened East 모델을 계산하는 C99 프로그램 모음입니다.
외부 라이브러리는 쓰지 않습니다.

| 프로그램 | 내용 | 입력 방식 |
|---|---|---|
| `east_ed` | 연속시간 tilted generator 정확 대각화: θ, k, e, χ, s\*(g, ε), strict 공존선 | 명령행 옵션 |
| `east_tps` | 이산시간(MCS) 동역학의 TPS (forward/backward shooting) | 입력 파일 |
| `east_tm` | 같은 MCS 동역학의 tilted 전이행렬 정확해: 유한 t_obs 모멘트, 장시간 극한 | 입력 파일 |
| `east_clone` | 이산시간 cloning (일정 개체 수) | 입력 파일 |
| `east_mbar` | 여러 TPS 창을 MBAR로 결합 | 입력 파일 |
| `east_reus` | K 창 replica-exchange umbrella sampling (REUS) | 입력 파일 |
| `east_edan` | `east_ed` peak 표로 가설 1–3 분석 (공존선 이동, L 스케일링, 두 준위 Δk, Δe) | 명령행 옵션 |
| `east_test` | 결정론적 검증 (824개) | — |
| `east_validate` | 시뮬레이션 대 정확해 통계 검증 (89개) | — |

## 1. 빌드와 검증

필요한 것: C99 컴파일러(gcc 또는 clang)와 make.

```bash
make              # bin/ 아래에 모든 프로그램 생성
make test         # 결정론적 검증, 약 15초
make validate     # 통계 검증, 약 2분
make OMP=1        # (선택) east_reus의 창들을 OpenMP로 병렬 실행
```

- Linux / macOS / WSL: 위 명령 그대로 사용
- Windows MinGW: `mingw32-make` 사용, 실행 파일은 `.\bin\east_tps.exe` 형식
- 모든 실행 파일은 프로젝트 최상위 폴더에서 `./bin/...`로 실행합니다

두 검증 모두 마지막 줄이 `0 failed`이면 정상입니다.

**`east_test` (결정론적, 824개)**

1–9. 연속시간 ED: 평형값, ε = 1 닫힌 식, k = 0 섹터 대 전체 공간, Lanczos 대 조밀 대각화, 유한차분, 맥스웰 관계, strict 공존선 닫힌 식, 곱 상태 정확해
10. 이산시간 전이행렬
   - L = 3, t_obs = 2에서 모든 시도 순서(격자점 선택과 뒤집힘 여부)를 직접 열거한 값과 ln Z, ⟨K⟩, ⟨E⟩가 10⁻¹²까지 일치
   - 평형값, 유한차분 대 정확 도함수, 유한 시간 → 장시간 극한
   - ε = 1, g = 0에서 θ_d = L ln(1 + μ(s)) 닫힌 식
11. REUS 시작 단계: `n_approach`를 일부러 작게 해서 못 들어간 창이 이웃 궤적 복사로 들어가는지, 모든 샘플과 마지막 궤적이 자기 창 안에 있는지 확인
12. REUS 체크포인트: 두 번 끊었다 이은 실행이 한 번에 끝까지 돌린 실행과 샘플, 카운터, 마지막 궤적까지 똑같은지, 매개변수가 다른 체크포인트를 거부하는지 확인

**`east_validate` (통계, 고정 시드, |z| ≤ 4이면 통과, 89개)**

| 항목 | 비교 대상 |
|---|---|
| 1 MCS 커널 | 고정된 x₀에서 2×10⁶번 진행한 x₁ 분포 대 정확한 M^L 행(χ² 검정), 평균 뒤집힘 수 |
| TPS | L = 8, t_obs = 40, (ε, s, g) 4가지: ⟨K⟩, ⟨E⟩, Var K, Var E, Cov 대 유한 시간 정확해 |
| Cloning | L = 8, 4가지: θ_d, k_d, e_d 대 전이행렬 최대 고유값 |
| MBAR | s = 0–0.04 다섯 창을 결합해 창 사이 s의 ⟨K⟩, ⟨E⟩, Var K 대 정확해 |
| Umbrella | L = 8, t_obs = 100, ε = 0.003, s = 0.03에서 K 창 8개 + MBAR: P(K < K_c), s = 0–0.06의 ⟨K⟩, Var K, 공존점(고정·자동 경계), 공존점에서 봉우리별 ⟨K⟩, ⟨E⟩ 대 정확한 P(K) |
| REUS | 같은 조건을 `east_reus`로: 위 항목 전부 + 교환 수락률, 왕복 횟수, 겹친 구간 KS 거리 |

## 2. 공통 규약

| 항목 | 정의 |
|---|---|
| 궤적 가중치 | e^{−sK−gE} (s > 0은 비활성, g > 0은 저에너지 궤적 선호) |
| 제약 | Cᵢ = n_{i−1} + ε(1 − n_{i−1}), 주기경계, 촉진 이웃은 왼쪽 |
| 평형 밀도 | c = 1/(1 + e^{1/T}), J = 1 |
| K | 전체 뒤집힘 수 (0→1, 1→0 모두) |

**연속시간 (`east_ed`)**: 전이율 w(0→1) = c·Cᵢ, w(1→0) = (1−c)·Cᵢ, E = ∫₀ᵗ Σᵢ nᵢ dt′.
θ는 계 전체·단위 시간당 값, k = −(1/L)∂θ/∂s, e = −(1/L)∂θ/∂g, χ_XY = (1/L)∂²θ/∂X∂Y.

**이산시간 (`east_tps`, `east_tm`, `east_clone`, `east_mbar`)**
- 1 MCS = L번 시도. 시도마다 격자점 i를 무작위로 고르고(복원추출), 확률 Cᵢ·c(nᵢ = 0일 때) 또는 Cᵢ·(1−c)(nᵢ = 1일 때)로 뒤집습니다. 이전 East TPS 코드의 갱신 규칙에 ε을 더한 것입니다.
- t_obs = MCS 수. 궤적은 슬라이스 x₀ … x_{t_obs}입니다.
- E = Σ_{t=0}^{t_obs} n(x_t) (t_obs + 1개 슬라이스 모두), K = t_obs번의 MCS 동안 일어난 뒤집힘 수
- 출력 정규화: k = ⟨K⟩/(L t_obs), e = ⟨E⟩/(L (t_obs+1)), χ_XY = Cov(X, Y)/(L t_obs)
- 장시간 극한 θ_d = ln λ_max(D_g^{1/2} M_s^L D_g^{1/2})는 MCS당·계 전체 값입니다.

이전 East TPS 코드와 다른 점은 셋입니다.
- ITERATION이 슬라이스 수(MCS는 ITERATION − 1회)였던 것을 t_obs = MCS 수로 바꿨습니다.
- `rand()` 대신 xoshiro256\*\* 난수를 씁니다. Windows의 `rand()`는 RAND_MAX가 32767이라 ε·c ≈ 3×10⁻⁴ 같은 작은 확률을 다루기에 해상도가 부족합니다.
- 입력 파일을 위치 순서 대신 키워드 형식으로 바꿨습니다.

## 3. 입력 파일 형식

```
# 주석
L      100        # 줄 끝 주석도 됩니다
eps    0.01
s      0.0:0.005:0.04   # 리스트를 받는 키워드: a,b,c 또는 시작:간격:끝
T      4/3              # 분수 가능
```

- 실행 파일 뒤에 경로를 주지 않으면 현재 폴더의 `system.inp`를 읽습니다: `./bin/east_tps`, `./bin/east_tps inputs/tps.inp`
- 프로그램이 쓰지 않는 키워드는 경고로 알려 줍니다 (오타 확인용).
- 예시 입력은 `inputs/` 폴더에 있습니다.

## 4. `east_tps` — TPS

이동은 shooting 하나뿐입니다. 확률 1/2로 forward(τ ∈ [0, t_obs−1]에서 그 뒤를 새 잡음으로 재생성) 또는 backward(τ ∈ [1, t_obs]에서 그 앞을 재생성)를 고르고, min(1, e^{−sΔK−gΔE})로 수락합니다. MCS 커널이 평형 곱 측도에 대해 가역이므로 backward 구간도 같은 커널로 만듭니다.

| 키워드 | 기본값 | 설명 |
|---|---|---|
| `L`, `eps`, `s`, `g`, `tobs`, `n_prod` | 필수 | |
| `T` | 1 | |
| `seed` | 0 | 0이면 시계 기반. 사용된 시드는 출력에 기록 |
| `n_equil` | 10000 | 샘플링 전 shooting 횟수 |
| `sample_every` | 10 | 몇 번 이동마다 (K, E)를 기록할지 |
| `init` | equilibrium | `equilibrium`(x₀를 평형에서 뽑고 앞으로 진행), `empty`(모든 슬라이스가 빈 비활성 궤적, ε > 0에서만 의미 있음), 또는 저장된 궤적 파일 |
| `save` | 없음 | 마지막 궤적 저장 (다음 창의 `init`으로 사용) |
| `output` | tps_samples.dat | 샘플 파일 |
| `nblocks` | 20 | jackknife 오차용 블록 수 |
| `progress` | 0 | 이 횟수마다 K, E를 stderr로 출력 |
| `umb_kmin`, `umb_kmax` | 없음 | umbrella 창: K를 [⌈umb_kmin·L·t_obs⌉, ⌊umb_kmax·L·t_obs⌋]로 제한 (k 단위로 입력) |
| `n_approach` | 1000000 | 시작 궤적이 창 밖일 때 창에 들어가기까지 허용하는 최대 이동 수 |

- umbrella 창이 있으면, 창 밖으로 나가는 이동은 기각합니다. 시작 궤적이 창 밖이면 먼저 K가 창에 가까워지는 이동만 받아 창 안으로 들어간 뒤(approach, 샘플에 포함 안 됨) `n_equil`, `n_prod`를 진행합니다. 이때 `.log`의 평균은 창 안으로 제한된 값이라 그 자체로는 물리량이 아니고, `east_mbar`로 합쳐야 합니다.
- 출력 파일: `# param` 헤더 줄들 → `move K E` 열 → 마지막에 `# summary` 줄들(평균, 분산, k, e, χ, 수락률, 블록 jackknife 오차). `east_mbar`가 이 헤더를 읽습니다.
- 궤적 파일: 슬라이스마다 `k_t 0101…` 한 줄 (텍스트, 플랫폼 무관)
- 속도(2코어 클라우드 기준, 이동 1회당 ms):

  | | t_obs = 200 | 400 | 600 |
  |---|---|---|---|
  | L = 25 | 0.015 | 0.028 | 0.047 |
  | L = 50 | 0.028 | 0.055 | 0.085 |
  | L = 100 | 0.058 | 0.114 | 0.165 |

여러 s 창을 이어서 돌리고 MBAR까지 하는 예시는 `bash scripts/run_tps_windows.sh`입니다. 이전 코드의 `run.sh`처럼 창마다 `system.inp`를 쓰고, 앞 창의 마지막 궤적을 다음 창의 초기값으로 씁니다.

### L, t_obs 스케일링: `scripts/run_tps_scaling.sh`

L × t_obs (× ε × g)의 모든 조합에 대해 s 창들을 이어서 돌리고 조합별로 MBAR까지 합니다.

1. 스크립트 위쪽의 변수를 정합니다: `EPS_LIST`, `G_LIST`, `L_LIST`, `TOBS_LIST`, `S_LIST`(창), `S_TARGET`(MBAR 격자), `N_PROD`, `INIT_FIRST`(첫 창의 시작 궤적: `equilibrium`, `empty`, 또는 궤적 파일), `OUT_ROOT`(결과 최상위 폴더, 기본 `tps`) 등. **`S_LIST`의 기본값은 예시일 뿐입니다.**
2. 실행 (프로젝트 최상위 폴더에서):

   ```bash
   bash scripts/run_tps_scaling.sh                    # 목록의 모든 조합을 차례로
   bash scripts/run_tps_scaling.sh 0.01 0.0 100 400   # 한 조합만: eps g L tobs
   ```

   조합마다 자기 폴더에만 쓰므로, 한 조합씩 따로 실행하면 병렬로 돌릴 수 있습니다 (터미널 여러 개, 또는 클러스터 작업 하나에 한 조합).
3. 결과는 `<OUT_ROOT>/eps<ε>_g<g>/L<L>_t<t_obs>/`에 생깁니다.
   - `s+0.0200.inp/.dat/.log/.traj`: 창별 입력, 샘플, 요약, 마지막 궤적
   - `mbar.inp`, `mbar.dat`: 이 조합의 MBAR 입력과 결과 (k, e, χ와 bootstrap 오차)

**계산 시간 추정**: 창 1개당 10⁶회 이동이면 9개 조합(L 3개 × t_obs 3개) 합계 약 10분입니다. 창 10개면 (ε, g) 하나당 약 1.7시간이고, (ε, g) 9개면 순차 실행으로 약 15시간입니다. 보고서의 "계 크기와 관측시간은 따로 스케일링"을 십자형(한 t_obs에서 L 3개 + 한 L에서 t_obs 3개 = 5개 조합)으로 하면 더 줄어듭니다.

**창(`S_LIST`) 고르는 법**
- 창은 공존점 s\* 주변에 둡니다. 연속시간 ED의 s\*는 시간 단위가 달라 참고값으로만 쓰고, 이산시간 위치는 예비 계산으로 확인하세요.
- 이웃 창의 K 분포가 겹쳐야 MBAR가 됩니다. 기준은 대략 Δs × √(Var K) ≲ 1–2이고, √(Var K)는 각 창 `.log`의 `Var K`에서 얻습니다. L·t_obs가 클수록 √(Var K)가 커지므로 창 간격을 좁혀야 합니다.
- 공존점 근처에서는 K 분포가 두 봉우리가 되어 Var K가 크게 늘어나므로, 거기서는 창을 더 촘촘히 둡니다.
- 큰 조합(예: L = 100, t_obs = 400) 하나를 `N_PROD`를 작게 해서 넓은 s 범위로 먼저 돌려 보고 창을 정하는 것을 권합니다.

**이력 현상(hysteresis) 점검**: 창을 큰 s부터 내려가며 돌려서 오름차순 결과와 비교합니다. 다음 세 가지를 바꾸세요.
- `S_LIST="0.040 0.035 0.030 0.025 0.020 0.015 0.010"` (내림차순)
- `INIT_FIRST=empty`: 첫 창(가장 큰 s)을 비활성 궤적에서 시작합니다. `equilibrium`으로 두면 내림차순도 활동적인 궤적에서 출발하므로 점검이 되지 않습니다.
- `OUT_ROOT=tps_down`: 오름차순 결과(`tps/`)를 덮어쓰지 않도록 다른 폴더에 씁니다.

같은 창의 k가 두 방향에서 오차 범위 안에서 같으면 출발 궤적의 영향이 없는 것이고, 다르면 그 창은 수렴하지 않은 것입니다.

**돌린 뒤 확인할 것**
- `.log`의 수락률: 너무 낮으면(수 % 이하) 상관이 길어 오차가 커집니다.
- `mbar.dat` 헤더의 `<K>_raw`와 `<K>_mbar`가 창마다 비슷해야 합니다.
- `ESS` 열이 작은 s(창 사이가 비어 있는 곳)의 값은 믿지 마세요.
- 공존점 판정에 쓸 P(k) 분포는 `mbar.inp`에 `hist_s`를 추가하고 `./bin/east_mbar <조합 폴더>/mbar.inp`를 다시 실행해서 봅니다.

### K 창 umbrella sampling: `scripts/run_umbrella.sh`

L·t_obs가 크면 활동/비활성 궤적군 사이 장벽 때문에 s 창 TPS가 한쪽에 갇힙니다(이력 현상). umbrella sampling은 s를 하나로 고정하고 K 범위를 겹치는 창들로 나눠 장벽 영역까지 직접 샘플링한 뒤, MBAR로 합쳐 P(K) 전체를 얻습니다. P_s'(K) ∝ P_s(K)·e^{−(s'−s)K}이므로 다른 s'로의 재가중이 정확하고, 공존점도 s' 방향 탐색으로 구합니다.

1. 스크립트 위쪽 변수를 정합니다.
   - `EPS`, `G`, `L`, `TOBS`
   - `S`: umbrella를 돌릴 s 하나. 공존점 근처로 둡니다
   - `KMIN_ALL`, `KMAX_ALL`: 창들이 덮을 k 범위. **두 봉우리를 모두 포함해야 합니다**
   - `NWIN`: 창 수. 이웃 창은 반씩 겹칩니다
   - `S_TARGET`: 재가중 s 격자이자 공존점 탐색 범위
   - `KCUT`: `auto`, k 값, 또는 `none`
   - `HIST_BINS`, `SMOOTH_RES`: 공존점 자동 경계의 구간 수 n과 평활 해상도 m
   - `HIST_S`: 재가중 P(k)를 출력할 s
   - `N_PROD`, `N_BOOT` 등
2. `bash scripts/run_umbrella.sh list`로 창 범위를 확인합니다.
3. 실행합니다 (프로젝트 최상위 폴더에서).

   ```bash
   bash scripts/run_umbrella.sh               # 모든 창을 차례로, 끝나면 MBAR
   bash scripts/run_umbrella.sh window 3      # 창 3만 (0부터 NWIN-1)
   bash scripts/run_umbrella.sh mbar          # 끝난 창들로 MBAR만
   ```

   창마다 `INIT`(기본 `empty`)에서 따로 출발하므로 창끼리 독립입니다. `window <i>`로 나눠 병렬 실행한 뒤 마지막에 `mbar`를 한 번 돌리면 됩니다.
4. 결과는 `<OUT_ROOT>/eps<ε>_g<g>/L<L>_t<t_obs>_s<s>/`에 생깁니다.
   - `w00.inp/.dat/.log/.traj`: 창별 입력, 샘플, 요약, 마지막 궤적
   - `mbar.dat`: s' 격자의 k, e, χ와 bootstrap 오차, 공존점(`# s_equal …` 줄), 봉우리별 평균(`# peaks …` 줄)
   - `mbar_hist.dat`: `HIST_S`의 재가중 P(k)

**창 설정 기준**
- 창 하나의 K 폭(= 2 × (KMAX_ALL − KMIN_ALL)/(NWIN+1) × L × t_obs)은 수백 정도가 적당합니다. 예비 계산(L = 50, t_obs = 600, 창 16개, 폭 약 210)에서 창 안 자기상관 시간이 샘플 10–150개로, 일반 TPS(수천)보다 훨씬 짧았습니다.
- 기본값(k 0–0.08, 창 22개)은 L = 100, t_obs = 600 기준이며 창 폭이 K 약 420입니다. L·t_obs가 작은 조합은 창 수를 줄여도 됩니다.

**점검할 것**
- `mbar.dat` 헤더에서 창마다 `<K>_raw`와 `<K>_mbar`(창 안 평균)가 비슷한지
- P(k)가 `KMIN_ALL`–`KMAX_ALL` 양 끝에서 충분히 작아지는지. 끝에서 크면 범위를 넓힙니다
- 공존점의 bootstrap 오차와, 성공한 replica 수(`n/N replicas`)

**계산 시간**: 창 1개 비용은 같은 L, t_obs의 일반 TPS와 같습니다. L = 100, t_obs = 600, 창 20개 × 10⁶회면 약 1시간입니다 (병렬이면 그만큼 단축).

### replica-exchange umbrella sampling: `east_reus`, `scripts/run_reus.sh`

K 창마다 따로 돌리는 umbrella sampling은 L·t_obs가 크면 이웃 창끼리 서로 맞지 않을 수 있습니다. 같은 K라도 궤적의 시간 구조(비활성 구간의 위치와 개수)가 다르고, 이 구조가 창 안에서 잘 바뀌지 않기 때문입니다. `east_reus`는 모든 창을 한 프로그램에서 동시에 돌리면서, `swap_every`회 이동마다 이웃 창의 궤적을 맞바꿉니다.

- 교환 조건은 두 궤적이 각각 상대 창 범위 안에 있는 것입니다. s, g가 같으므로 그러면 항상 수락됩니다.
- 궤적이 창 사이를 오가며 구조가 쉽게 바뀌는 쪽(아주 비활성이거나 아주 활동적인 쪽)에서 바뀐 뒤 돌아올 수 있습니다.
- 출력 창 파일 형식이 `east_tps`와 같아서 `east_mbar`를 그대로 씁니다. 네 번째 열은 그 샘플을 만든 replica 번호입니다.

| 키워드 | 기본값 | 설명 |
|---|---|---|
| `L`, `eps`, `s`, `g`, `tobs`, `nwin`, `k_lo`, `k_hi`, `n_prod` | 필수 | 창 i = [k_lo + iD, k_lo + (i+2)D], D = (k_hi − k_lo)/(nwin+1) |
| `T`, `seed` | 1, 0 | |
| `init` | empty | `empty`, `equilibrium`, 또는 `%02d`가 들어간 궤적 파일 패턴(창 번호). 없는 파일은 빈 궤적으로 시작 |
| `n_approach` | 1000000 | 각 replica가 창에 들어가기까지 한 번 시도할 때의 최대 이동 수 (아래 시작 단계 참고) |
| `n_equil` | 20000 | 창당 평형화 이동 수 (교환 포함) |
| `swap_every`, `sample_every` | 10, 10 | 교환 간격, 샘플 간격 (이동 수, sample_every는 swap_every의 배수) |
| `progress` | 0 | 이 주기(cycle)마다 진행 상황 출력 |
| `output_prefix` | reus_ | 출력 파일 앞부분 (폴더 포함 가능) |
| `save_traj` | 1 | 창별 마지막 궤적 저장 여부 |
| `checkpoint` | 없음 | 체크포인트 파일. 주면 중간 저장과 이어서 계산을 함 (아래 참고) |
| `checkpoint_every` | 1800 | 체크포인트 간격 (실제 경과 시간, 초) |
| `stop_after` | 0 | 이 시간(초)이 지나면 체크포인트를 쓰고 멈춤 (0: 안 멈춤) |

**시작 단계 (approach)**
- 각 replica는 `init` 궤적에서 출발해, 교환 없이 혼자 자기 창에 들어갑니다. 창까지 K 거리가 줄어드는 이동만 받습니다.
- 이 방식은 창을 건너뛰어 반대편으로 갈 수 있습니다. 특히 빈 궤적에서 한 번에 활동 상태로 넘어가면, 두 봉우리 사이 골 근처 창에서는 K를 다시 줄이는 이동이 드물어 오래 걸리거나 `n_approach` 안에 못 들어갈 수 있습니다.
- 못 들어간 창은, 이미 들어간 창 가운데 K가 그 창에 가장 가까운 replica의 궤적을 복사해 다시 시작합니다. 이웃 창끼리 절반씩 겹치므로 복사한 궤적은 대개 이미 창 안이거나 아주 가깝습니다.
  - 이 과정은 모든 창이 들어갈 때까지 반복합니다. 연속 세 번 새로 들어간 창이 없으면 멈추고 `error: window … not reached`를 냅니다.
  - 복사가 일어나면 화면에 `note: window … restarting from the trajectory of window …`가 나옵니다.
  - 복사한 replica는 처음에 이웃과 같은 궤적이지만 난수 흐름이 달라 곧 갈라지고, 이어지는 `n_equil` 동안 교환과 함께 평형화됩니다.

**체크포인트와 이어서 계산**
- `checkpoint`를 주면 전체 상태(궤적, 난수 상태, 창 배치, 카운터, 지금까지의 샘플)를 저장합니다.
  - 저장 시점: 시작 단계가 끝난 직후, `checkpoint_every`초마다, 멈추라는 요청을 받았을 때
  - 멈추라는 요청: `stop_after`초 경과, 또는 SIGTERM, SIGUSR1 신호
  - `<파일>.tmp`에 먼저 쓰고 이름을 바꾸므로, 쓰는 도중에 죽어도 이전 체크포인트는 남습니다.
- 멈추면 종료 상태 3으로 끝납니다. 같은 입력으로 다시 실행하면 체크포인트에서 이어서 계산합니다.
- 중간에 여러 번 끊었다 이어도 결과는 한 번에 끝까지 돌린 것과 비트 단위로 같습니다(`east_test` 12번과 L = 100 실제 조건으로 확인).
- 매개변수가 다른 체크포인트는 거부합니다. 처음부터 다시 하려면 체크포인트 파일을 지우세요.
- 모든 출력 파일을 쓴 뒤 체크포인트 파일을 지웁니다.
- 시작 단계 도중에 멈추라는 요청을 받으면 체크포인트 없이 멈춥니다. 다음 실행은 처음부터 시작합니다.

출력 파일은 다음과 같습니다.
- `<prefix>w00.dat …`: 창별 샘플 (`east_mbar` 입력)
- `<prefix>w00.traj …`: 창별 마지막 궤적. 다음 실행의 `init` 패턴으로 다시 쓸 수 있습니다
- `<prefix>summary.dat`(같은 내용이 화면에도 나옴):
  - 창별 ⟨K⟩, 수락률, 독립 샘플 수(N_eff), 창 진입에 걸린 이동 수(`approach`, 모든 시도 합), 복사해 온 창 번호(`from`, 복사 안 했으면 −1)
  - 이웃 쌍별 교환 수락률과 겹친 구간 KS 거리
  - replica 왕복 횟수(창 0 → 맨 위 → 창 0)

**실행** (프로젝트 최상위 폴더에서; 서버에서는 아래 SLURM 제출을 쓰세요):

```bash
make OMP=1                                     # 병렬 실행용 빌드 (한 번)
OMP_NUM_THREADS=4 bash scripts/run_reus.sh     # REUS(또는 이어서) 후 MBAR, 공존점
EPS=0.003 G=0.05 bash scripts/run_reus.sh      # 매개변수를 환경 변수로 바꿔서 실행
bash scripts/run_reus.sh mbar                  # MBAR만 다시
bash scripts/run_reus.sh params                # 실제로 쓸 매개변수 출력
```

- `scripts/run_reus.sh` 위쪽의 모든 매개변수는 같은 이름의 환경 변수가 있으면 그 값을 씁니다.
- 체크포인트는 항상 `<결과 폴더>/reus.ckpt`에 둡니다(`CHECKPOINT_EVERY`, `STOP_AFTER`로 조절). 멈춘 뒤 같은 명령을 다시 실행하면 이어서 계산합니다. 이때 종료 상태는 3입니다.
- REUS가 끝나면 그 매개변수를 `<결과 폴더>/reus.done`에 남깁니다.
  - 같은 매개변수로 다시 실행하면 REUS를 건너뛰고 MBAR만 다시 합니다.
  - 매개변수가 다르면 에러로 멈춥니다. `OUT_ROOT`를 바꾸거나 `reus.done`을 지우세요.

결과는 `<OUT_ROOT>/eps<ε>_g<g>/L<L>_t<t_obs>_s<s>/`에 생깁니다. OpenMP 빌드는 직렬 빌드와 결과가 비트 단위로 같고(검증함), 2코어에서 약 1.7배 빨랐습니다. 전체 계산량은 창별 umbrella와 같으므로, L = 100, t_obs = 600, 창 22개 × 10⁶회면 직렬로 약 1시간이고 코어 수만큼 줄어듭니다.

**기본값**은 L = 100, t_obs = 600, ε = 0.01, g = 0 기준입니다: k 0–0.11, 창 30개(창 폭 k 약 0.007), s = 0.018. 예비 REUS(창당 10⁵회)에서 공존점 근처(s ≈ 0.0175)의 활동 봉우리가 k ≈ 0.07–0.08에 있어 k 0.08에서 잘렸기 때문에, 위쪽 끝을 0.11로 넓혔습니다.

**수렴 판단**
- 모든 이웃 쌍의 교환 수락률이 0보다 충분히 커야 합니다(수 % 이상). 0인 쌍이 있으면 그 사이가 끊긴 것이니 창을 더 촘촘히 둡니다.
- 왕복 횟수가 replica마다 여러 번 이상이어야 창 전체가 섞인 것입니다.
- 겹친 구간 KS 거리가 작아야 합니다. 검증 계산에서는 0.012 이하였습니다.

### SLURM 제출: `scripts/submit_reus.sh`

서버에서는 REUS를 SLURM 작업으로 제출합니다. 제출 스크립트가 아래 조건을 맞춰 줍니다.

| 설정 (`scripts/submit_reus.sh` 위쪽) | 기본값 | 의미 |
|---|---|---|
| `NODE` | node3 | 계산노드 3의 SLURM 이름. 서버 구성이 바뀌면 `sinfo -N -l`로 다시 확인합니다 |
| `PARTITION` | 비어 있음 | 파티션 (비우면 클러스터 기본값) |
| `CPUS` | 4 | 작업당 코어 수 = OpenMP 스레드 수 |
| `TIME` | 24:00:00 | 한 번 실행에서 요청할 최대 시간. 파티션 최대 시간이 더 짧으면 자동으로 그 시간을 사용 |
| `JOB_NAME` | soft_east | 이 이름의 작업은 `--dependency=singleton`으로 한 번에 하나씩만 돕니다. 그래서 동시에 쓰는 코어는 최대 `CPUS`개입니다 |
| `MAX_REQUEUE` | 10 | 작업 하나가 자동으로 requeue될 수 있는 최대 횟수 |
| `MARGIN` | 1800 | 시간 제한보다 이만큼(초) 먼저 REUS가 체크포인트를 쓰고 멈춤 |
| `SIGNAL_LEAD` | 600 | 시간 제한보다 이만큼(초) 먼저 SLURM이 작업 스크립트에 SIGUSR1을 보냄 |

```bash
bash scripts/submit_reus.sh                              # run_reus.sh 기본값으로 제출
bash scripts/submit_reus.sh EPS=0.003 G=0.05 TOBS=400    # 이 작업만 매개변수 변경
bash scripts/submit_reus.sh MODE=mbar EPS=0.003 G=0.05   # MBAR만 (REUS가 끝난 결과에)
bash scripts/submit_reus.sh --dry-run EPS=0.003          # 제출하지 않고 내용만 확인
bash scripts/submit_reus.sh --from slurm_jobs/<파일>.env  # 예전 작업과 같은 매개변수로 다시 제출
```

- **매개변수 고정:** 제출할 때 모든 매개변수를 `slurm_jobs/<이름>.env`에 저장하고 작업은 그 값을 씁니다. 그래서 제출 뒤에 `run_reus.sh`를 고쳐도 대기 중인 작업은 바뀌지 않습니다. 이름이 틀린 매개변수는 거부합니다.
- **시간 제한 자동 조정:** 제출할 때 `sinfo`로 지정한 파티션(비어 있으면 기본 파티션)의 최대 시간을 읽습니다. `TIME`보다 짧으면 파티션 최대 시간을 `--time`에 쓰고, 체크포인트 중단 시점과 SIGUSR1 사전 통보 시간도 그 안에 맞게 줄입니다. 기본 파티션을 찾지 못하면 `PARTITION`을 지정해야 합니다.
- **제출 목록:** 제출할 때마다 `slurm_jobs/jobs.txt`에 `작업번호 제출시각 모드 결과폴더 매개변수파일`을 한 줄씩 남깁니다.
- **로그:** `slurm_logs/soft_east_<작업번호>.out`에 쌓입니다. requeue된 실행도 같은 번호, 같은 파일에 이어서 씁니다.
- **24시간 안에 안 끝날 때 (자동 requeue)**
  1. REUS는 시작 후 `TIME − MARGIN`(기본 23시간 30분)이 지나면 체크포인트를 쓰고 멈춥니다.
  2. 작업 스크립트(`scripts/slurm_reus.sh`)가 `scontrol requeue`로 같은 작업을 다시 대기열에 넣습니다.
  3. 다시 시작하면 체크포인트에서 이어서 계산하고, 끝나면 MBAR까지 합니다.
  - 그때도 다른 일(예: MBAR)을 하고 있으면, 시간 제한 10분 전 SIGUSR1을 받은 작업 스크립트가 REUS에 체크포인트를 쓰게 한 뒤 requeue합니다. REUS가 끝난 뒤라면 다음 실행은 MBAR만 다시 합니다.
  - requeue가 `MAX_REQUEUE`번을 넘으면 더 하지 않고 로그에 `not requeued`를 남깁니다. 이어서 하려면 `--from`으로 같은 매개변수를 다시 제출하면 됩니다. 체크포인트가 있으면 이어서 계산합니다.
- **클러스터 조건:** 사용자가 자기 작업을 `scontrol requeue`할 수 있어야 합니다. 처음 제출한 작업의 로그에서 requeue 줄을 확인하세요.

**작업 확인과 취소**

```bash
squeue -u $USER -n soft_east                          # 대기/실행 중인 작업
tail -n 30 slurm_logs/soft_east_<작업번호>.out         # 진행 상황 (cycle …, requeue …)
sacct -j <작업번호> --format=JobID,State,Elapsed,ExitCode,NodeList
scontrol show job <작업번호>                        # sacct가 비활성화된 서버의 대체 확인
cat slurm_jobs/jobs.txt                               # 제출 목록
scancel <작업번호>                                     # 취소 (체크포인트는 남아 다음 제출 때 이어짐)
```

이 서버는 현재 SLURM accounting storage가 비활성화되어 `sacct`가 실패합니다. 작업 종료 직후에는 `scontrol show job <작업번호>`의 `JobState`, `ExitCode`, `Restarts`를 확인합니다.

## 5. `east_tm` — 이산시간 정확해 (작은 L)

| 키워드 | 설명 |
|---|---|
| `mode` | `finite`: t_obs MCS 궤적 앙상블의 정확한 ln Z, ⟨K⟩, ⟨E⟩, Var, Cov (TPS 검증용). `asym`: 장시간 극한 θ_d, k_d, e_d, λ₂/λ₁, χ_KK (cloning 검증, 연속시간과 경향 비교용). `pk`: 정확한 P(K) 분포와, `kcut`을 주면 공존점 (umbrella 검증용) |
| `L`, `eps`, `s`, `g` | 필수, `s`와 `g`는 리스트 가능 |
| `tobs` | `finite`, `pk`에서 필수 |
| `kcut`, `hist_bins`, `smooth_res` | `pk`에서 공존점 계산 (`east_mbar`와 같은 정의) |
| `basis` | `k0`(기본) 또는 `full` |
| `fd_h` | `asym`의 χ_KK = −dk/ds 중심차분 간격 (기본 1e-4, 0이면 생략) |
| `output` | 기본 stdout |

- L ≤ 14 정도가 실용 범위입니다 (L = 12 `finite`, t_obs = 600이 약 0.1초).
- 도함수는 순방향 자동미분으로 정확히 계산합니다 (유한차분 아님).
- `pk`는 계산량이 (L·t_obs)²에 비례해서 L·t_obs가 수천 이하일 때만 실용적입니다 (L = 8, t_obs = 200이 약 0.6초).
- `asym`에서 ε = 0이면 빈 격자 상태를 뺀 섹터를 쓰고, 빈 상태는 θ = 0입니다.

## 6. `east_clone` — cloning

이전 cloning 코드와 같은 이산시간 방식입니다. 각 clone이 MCS를 1회 진행하고 가중치 e^{−sk−gn}을 받으며, 개체 수를 일정하게 유지하는 systematic resampling을 합니다.

| 키워드 | 기본값 | 설명 |
|---|---|---|
| `L`, `eps`, `s`, `g`, `n_clones`, `tobs` | 필수 | `s`, `g`는 리스트 가능 |
| `t_skip` | 100 | 처음 이만큼의 단계는 모든 추정에서 제외 |
| `t_tail` | 100 | 마지막 이만큼의 단계는 k, e 이력에서 제외 |
| `n_runs` | 10 | 독립 반복 횟수 (오차 = 반복 간 표준오차) |
| `seed`, `print_runs`, `output` | | |

- θ_d = 단계별 ln⟨w⟩의 평균 (t_skip 이후)
- k, e는 마지막 개체들이 이어받은 이력에서 구합니다. 마지막 시각 근처는 편향 궤적의 내부 통계와 다르므로 `t_tail`만큼 뺍니다. **`t_tail`은 이완 시간보다 길어야 합니다.** 검증에서 L = 8, ε = 0.01 비활성 쪽은 t_tail = 0이면 k가 24% 크게 나왔고, t_tail = 300이면 정확해와 맞았습니다.

## 7. `east_mbar` — MBAR

상태 k의 편향을 u_k = s_k K + g_k E (umbrella 창 [Kmin_k, Kmax_k] 밖이면 +∞)로 두고, 볼록 목적함수에 Newton법을 써서 자유에너지를 구합니다. 초기값은 이웃 umbrella 창의 겹침 비율(또는 s 방향 열역학 적분)로 잡습니다.

| 키워드 | 설명 |
|---|---|
| `window <파일>` | `east_tps` 샘플 파일, 한 줄에 하나씩 반복. (s, g, umbrella 창)이 같은 파일은 한 상태로 합칩니다. s 창 TPS와 umbrella 창을 섞어도 됩니다 |
| `target_s`, `target_g` | 결과를 낼 (s, g) 격자 |
| `skip` | 각 파일 앞부분에서 버릴 샘플 수 (기본 0) |
| `n_boot`, `block`, `seed` | moving-block bootstrap 오차 (기본 0 = 오차 없음, block = 10 샘플) |
| `kcut` | 선택. 공존점(두 봉우리 무게가 같은 s)을 target_s 범위에서 찾음. 숫자면 P(k < kcut) = 1/2, `auto`면 아래 정의 |
| `hist_bins` | 히스토그램 구간 수 n (기본 60). 경계 위치의 해상도를 정함 |
| `smooth_res` | `auto`의 평활 해상도 m (기본 30). 각 구간을 앞뒤 round(n/(2m))개 구간과 평균함. 평균 창은 히스토그램 범위의 약 1/m(+ 구간 1개)이라 구간 수를 바꿔도 거의 일정함 (0이면 평활 안 함) |
| `hist_s`, `hist_output` | 선택. 재가중된 P(k) 분포 출력 |
| `tol`, `output` | |

- 출력 열: `g s k e chi_KK chi_EE chi_KE rho ESS` (+ bootstrap 시 오차 5열)
- 헤더에 상태별 자유에너지 f와, 각 창의 원래 ⟨K⟩와 MBAR로 다시 계산한 ⟨K⟩가 나옵니다 (일치하면 정상).
- **ESS**(Kish 유효 샘플 수)가 작은 s는 창들이 그 영역을 덮지 못한 것이니 믿지 마세요.
- 모든 창 파일의 L, t_obs, ε, T가 같아야 합니다.
- **`kcut auto`의 정의**
  - 재가중 P(k) 히스토그램(구간 n개)을 이동평균합니다. 각 구간을 앞뒤 round(n/(2m))개 구간과 평균하므로, 평균 창은 히스토그램 범위의 약 1/m이고 구간 수와 거의 무관합니다. 구간 수는 n개 그대로 유지됩니다.
  - 가장 높은 봉우리와, 둘 사이 골이 가장 깊은 다른 봉우리를 고릅니다. 골 깊이가 낮은 쪽 봉우리 높이의 10% 이상이어야 두 봉우리로 인정합니다.
  - 두 봉우리 사이 최저점을 경계로 무게가 반반인 s를 찾고, 그 s에서 경계를 다시 구하는 과정을 수렴할 때까지 반복합니다.
  - 공존점에서 분포가 두 봉우리가 아니면 `nan`이 나옵니다(작은 계에서 흔함).
  - 이 정의는 제가 정한 것이니, 보고서에 쓰기 전에 동의하시는지 확인해 주세요. 숫자 `kcut`과 결과를 비교해 보는 것도 방법입니다.
- **봉우리별 평균 (`# peaks` 줄, 가설 3용)**
  - `kcut`을 주면 각 target_g의 공존점 s\*에서, 경계 아래(비활동, k < kcut)와 위(활동, k ≥ kcut)로 나눈 조건부 평균을 냅니다.
  - 열: `g s* kcut k_in k_act dk e_in e_act de sdk_in sdk_act`
    - dk = k_act − k_in, de = e_act − e_in
    - sdk_in, sdk_act: 봉우리 안에서 k의 표준편차
    - e = ⟨E⟩/(L(t_obs+1))로, 위 표의 e와 같은 정의입니다.
  - bootstrap을 켜면 `err_k_in err_k_act err_dk err_e_in err_e_act err_de nrep`가 붙습니다.
    - 각 bootstrap 복제본에서 s\*와 경계를 다시 구한 뒤 그 점에서 계산하므로, 오차에 s\*와 경계의 불확실성이 들어갑니다.
    - nrep는 s\*를 찾은 복제본 수입니다.
  - 공존점이 `nan`이면 이 줄도 `nan`입니다.
  - 데이터 줄만 뽑기: `grep '^# peaks ' mbar.dat`
  - 해석할 때 주의할 점
    - 분포가 한 봉우리여도 경계로 나누면 dk > 0이 나옵니다. 정규분포를 중앙값에서 나누면 dk ≈ 1.6σ입니다. 그래서 dk는 sdk_in, sdk_act와 같이 봐야 합니다.
    - 유한 t_obs에서는 비활동 봉우리의 궤적도 시간 경계 근처가 활동적이라 k_in이 0보다 꽤 큽니다. 따라서 dk, de는 t_obs = 200, 400, 600 결과를 비교해 t_obs 의존성을 봐야 합니다.

## 8. `east_ed` — 연속시간 정확 대각화

```
./bin/east_ed <mode> [options]
```

| 모드 | 내용 |
|---|---|
| `scan` | (s, g) 격자 위의 모든 관측량 |
| `peak` | 각 g에서 χ_KK 최대점 s\* (ε > 0 필요) |
| `strict` | ε = 0에서 θ_nonempty(s\*, g) = 0인 s\* |

리스트 인자는 `0.1`, `8,10,12`, `-0.15:0.05:0.15` 형식을 받습니다.

| 옵션 | 기본값 | 설명 |
|---|---|---|
| `-L` | 10 | 계 크기 (실용 범위 ≤ 18) |
| `--T` | 1 | 온도 |
| `--eps` | 0.01 | softness |
| `--g` | 0 | 에너지 편향 |
| `--s` | -0.05:0.005:0.10 | s 값, `peak`에서는 거친 탐색 격자 |
| `--basis` | k0 | `k0` 또는 `full` |
| `--chi` | lr | `lr`(선형응답), `fd`(유한차분), `none` |
| `--fd-h` | 1e-4 | 유한차분 간격 |
| `--stol` | 1e-7 | s\* 허용오차 |
| `--side δ` | 없음 | `peak`에서 s\* ∓ δ의 k, e와 Δk, Δe 출력 |
| `--tol`, `--mmax` | 1e-11, 100 | Lanczos 설정 |
| `--out` | stdout | 출력 파일 |

```bash
./bin/east_ed scan -L 12 --eps 0.01 --g 0.05 --s 0.0
./bin/east_ed peak -L 8,10,12,14 --eps 0.003,0.01 --g -0.15:0.05:0.15 --s -0.5:0.005:0.8 --out results/peak.dat
./bin/east_ed strict -L 8,10,12,14 --g -0.15:0.05:0.15 --out results/strict.dat
```

보고서 계획의 ED 계산 전체는 `bash scripts/run_ed_plan.sh`로 돌릴 수 있습니다 (2코어에서 약 25분). 결과는 `results/`에 있습니다.

- 출력 열 `scan`: `L eps g s theta theta2 gap k e chi_KK chi_KE chi_EE rho resid nmv cg_it ok`
- 출력 열 `peak`: `L eps g s_star chi_KK_max theta gap k e chi_KE chi_EE rho edge ok`
- 출력 열 `strict`: `L g s_star s_exact k_active e_active slope gap_nonempty resid ok`
- 수치 방법: k = 0 병진 섹터, 상세 균형으로 대칭화, thick-restart Lanczos, Hellmann–Feynman 1차 도함수, 축소 레졸벤트(CG) 2차 도함수

### 가설 분석용 추가 계산: `scripts/run_ed_fine.sh`

`results/peak_T1.dat`(ε 5개, L = 8–14)만으로는 ε 간격이 넓고 L이 작아서 아래 계산을 더 합니다 (2코어에서 약 40분).

| 단계 | 내용 | 결과 |
|---|---|---|
| 1 | ε = 0.0005, 0.0015, 0.002, 0.004, 0.005, 0.007, 0.015, 0.02, L = 8–14, s 범위 `-0.5:0.005:0.8` | `results/peak_T1_fine.dat` |
| 2 | L = 16, 모든 (ε, g) | `results/peak_T1_L16.dat` |
| 3 | L = 18, ε = 0.0003, 0.001, 0.003, 0.01, 0.03 × g = −0.15, 0, 0.15 | `results/peak_T1_L18.dat` |
| 4 | strict 공존선, L = 16, 18 | `results/strict_T1_L16_18.dat` |

```bash
bash scripts/run_ed_fine.sh            # 1–4단계 전부
bash scripts/run_ed_fine.sh 2          # 한 단계만
NPROC=4 bash scripts/run_ed_fine.sh    # 동시에 돌릴 작업 수 (기본: 코어 수)
bash scripts/run_ed_fine.sh merge      # 작업별 파일(results/tmp_ed_fine/)로 표만 다시 만들기
```

- L ≥ 16은 s 전체 범위를 훑기에는 너무 느려서, 같은 (ε, g)의 L = 14 s\* 주변 [s\* − 0.02, s\* + 0.03] (간격 0.001)만 거칠게 탐색합니다. 그래서 L ≥ 16 행의 `edge` 판정은 이 창 기준입니다.
- L = 14에서 `nan`인 점(ε = 0.02, 0.03, g = −0.15)은 L ≥ 16에서 계산하지 않습니다.
- 한 점에 걸리는 시간: L = 16 약 40초, L = 18 약 2–3분.

### 가설 분석: `east_edan`, `scripts/run_ed_analysis.sh`

```bash
./bin/east_edan <mode> [옵션] <peak 표 ...>
bash scripts/run_ed_analysis.sh        # 아래 네 모드를 모두 실행 → results/analysis/
```

| 모드 | 계산 | 가설 | 결과 파일 |
|---|---|---|---|
| `h2` | 각 L의 r = (s\*_soft − s\*_strict)/ε, 가장 큰 L 세 개로 L → ∞ 외삽, 예측값 1/[2c(1−c)] | 2 | `h2_T1.dat` |
| `h1` | 연속한 L 쌍의 α = d ln χ_KK,max / d ln L, β = d ln(L·gap) / d ln L | 1 | `h1_T1.dat` |
| `cross` | g와 L 쌍마다 α가 기준값을 지나는 ε, β = 0이 되는 ε (ln ε에서 선형 보간) | 1 | `cross_T1.dat` |
| `h3` | 두 준위 근사 Δk₂ = √(2χ_KK·gap/L), Δe₂ = √(2χ_EE·gap/L), strict 활동 가지 k, e와의 비 | 3 | `h3_T1.dat` |

| 옵션 | 기본값 | 의미 |
|---|---|---|
| `--alpha a,b,..` | 0.25,0.5,1 | `cross`의 α 기준값 (스크립트는 0.25,0.5,1,1.75) |
| `--strict <file>` | 없음 | `h3`에 쓸 strict 표, 여러 번 줄 수 있음 |
| `--out <file>` | stdout | 출력 파일 |

- `edge = 1`, `ok = 0`, s\* = nan인 행은 없는 점으로 취급합니다. 같은 점이 여러 표에 있으면 먼저 준 표의 값을 씁니다.
- s\*_strict는 표 머리말의 c로 닫힌 식 e^{−2s\*} = [(1+2g)² − (1−2c)²]/[4c(1−c)]에서 계산합니다.
- `h2` 외삽
  - `r_geo`: 등간격 L 세 개의 차이가 등비로 줄어든다고 보고 외삽 (Aitken)
  - `r_pow`: r(L) = r_∞ − A L^{−p}를 세 점에 맞춤 (`p` 열)
  - 두 값의 차이를 외삽 불확실성의 대략적인 크기로 봅니다. 세 값이 단조롭지 않으면 `nan`입니다.
- `cross`의 `where` 열
  - `0`: 교차 1번, `2`: 여러 번 (가장 작은 ε)
  - `1`: 모든 ε에서 기준보다 큼, `-1`: 모든 ε에서 기준보다 작음 → 교차점이 계산한 ε 범위 밖
  - L = 18은 ε 5개뿐이라 16–18 쌍의 교차값은 거친 보간입니다.

**지표 읽는 법**

- 두 상태가 섞이는 구간(1차 전이 쪽)에서는 χ_KK,max ≈ L·Δk²/(2·gap)입니다.
  - gap이 L에 따라 닫히면 α > 1이고 L이 커질수록 커집니다. crossover면 χ가 포화해 α → 0입니다.
  - β < 0이면 gap이 1/L보다 빨리 닫히는 것이고, β → 1이면 gap이 L과 무관한 것입니다.
- 위 관계 때문에 α ≈ 2 − β + 2 d ln Δk / d ln L입니다. 즉 α와 β는 서로 독립인 지표가 아닙니다.
- 두 준위 근사(`h3`)는 rho = χ_KE/√(χ_KK χ_EE)가 1에 가까울 때만 의미가 있습니다. 현재 표에서 rho는 ε ≤ 0.01이면 0.977 이상, ε = 0.015–0.02이면 0.965 이상, ε = 0.03이면 0.95–0.97입니다.

**수정한 작업가설 (2026-11-17)**

1. `ε_c(g)`는 g가 커질수록 커지며, g = 0에서는 10⁻³ 정도이다. 시뮬레이션에서는 주어진 `(L, t_obs)`에서 χ_KK가 최대인 s의 P(k)가 두 봉우리인 가장 큰 ε를 유한계 `ε_c`로 정의하고, L과 t_obs 의존성을 함께 본다.
2. 작은 ε에서 `s*(ε,g) = s*_strict(g) + A(g) ε + O(ε²)`이다. `A(0) = 1/[2c(1−c)] = 2.543`이고 A(g)는 g가 커질수록 작아진다. 2.543은 g = 0, ε → 0 극한의 계수이지 모든 g와 ε에 적용하는 값이 아니다.
3. ε가 아래에서 `ε_c`로 갈수록 공존 봉우리의 Δk와 Δe가 연속적으로 줄어든다. 유한계에서는 두 봉우리가 유한한 폭을 가진 채 합쳐지므로 Δk = 0 자체보다 Δk를 `sdk_in`, `sdk_act`와 함께 비교해 분리도가 줄어드는지 검정한다.

**현재 결과 (T = 1, L ≤ 18, 2026-11-17)**

가설 2: `A(g) ≈ Δs*/ε`의 L → ∞ 외삽값 (`r_geo`–`r_pow`; g = 0 예측 1/[2c(1−c)] = 2.543)

| g | ε = 0.0003 | ε = 0.003 | ε = 0.01 | ε = 0.03 |
|---|---|---|---|---|
| −0.15 | 3.91–4.05 (L에 대해 아직 수렴 안 함) | 2.46 | 1.47 | nan |
| −0.10 | 3.29–3.40 | 2.52–2.54 | 1.71 | 0.81 |
| −0.05 | 2.85–2.90 | 2.44–2.47 | 1.78–1.79 | 1.00 |
| 0 | 2.54–2.55 | 2.33–2.35 | 1.80 | 1.10 |
| 0.05 | 2.31–2.33 | 2.19–2.23 | 1.78–1.79 | 1.15 |
| 0.10 | 2.15–2.16 | 2.08–2.11 | 1.75–1.76 | 1.18 |
| 0.15 | 2.02–2.03 | 1.99–2.01 | 1.72–1.73 | 1.20 |

- g = 0, ε → 0에서는 예측 2.543과 맞습니다.
- g ≠ 0이면 ε → 0에서도 계수가 g에 따라 달라집니다 (g = 0.15: 2.03, g = −0.1: 3.3–3.4).
- 작은 ε의 A(g)는 g가 커질수록 작아집니다 (g = −0.15에서 약 4.0, g = 0.15에서 약 2.0).
- ε가 커지면 계수가 줄어듭니다. g = 0, ε = 0.01에서 1.80으로, 선형 예측보다 약 30% 작습니다 (s\* ≈ 0.0180).
- 작업 해석: g = 0에서는 비활동 가지의 1차 보정과 평형 활동 가지로부터 `A(0) = c/k_eq = 1/[2c(1−c)]`가 나옵니다. g ≠ 0에서는 활동 가지가 평형에서 벗어나 그 가지에도 ε 보정이 생기므로 A(g)가 달라진다고 봅니다. 이는 현재 수치와 맞는 해석이며 별도 유도·검증이 필요합니다.

가설 1: L 쌍 14–16에서 기준을 지나는 ε (g = −0.15 … 0.15)

| 기준 | −0.15 | −0.10 | −0.05 | 0 | 0.05 | 0.10 | 0.15 |
|---|---|---|---|---|---|---|---|
| α = 0.25 | 6.1e-4 | 1.5e-3 | 2.6e-3 | 3.7e-3 | 4.8e-3 | 6.0e-3 | 7.0e-3 |
| α = 0.5 | 3.0e-4 | 8.9e-4 | 1.6e-3 | 2.5e-3 | 3.3e-3 | 4.1e-3 | 4.9e-3 |
| α = 1 | < 3e-4 | 3.4e-4 | 7.7e-4 | 1.3e-3 | 1.8e-3 | 2.4e-3 | 3.0e-3 |
| α = 1.75 | < 3e-4 | < 3e-4 | < 3e-4 | 4.0e-4 | 7.0e-4 | 1.1e-3 | 1.4e-3 |
| β = 0 | < 3e-4 | < 3e-4 | < 3e-4 | 3.6e-4 | 6.7e-4 | 1.0e-3 | 1.4e-3 |

- 어느 기준, 어느 L 쌍에서도 교차 ε는 g가 커질수록 커집니다.
- 하지만 값 자체는 기준에 따라 약 10배 차이가 납니다. 또 L에 따라 움직입니다: α = 0.25, 0.5는 대체로 L이 커지면 작아지고, α = 1.75와 β = 0은 커집니다.
- g ≥ 0에서는 α = 1 교차점만 L 쌍 10–12부터 16–18까지 거의 그대로입니다 (g = 0: 1.1–1.3e-3, g = 0.15: 2.8–3.1e-3).
- g = 0에서 α = 1.75, β = 0, α = 1 같은 강한 1차 전이 판정 기준은 약 0.4–1.3×10⁻³입니다. α = 0.25, 0.5 같은 느슨한 crossover 기준은 더 큰 값을 줍니다.
- 그러므로 ED(L ≤ 18)로는 ε_c의 정확한 값보다 g에 대한 경향과 10⁻³ 크기 범위를 말할 수 있습니다. 최종 판정은 위의 유한 `(L, t_obs)` REUS 정의로 합니다.

가설 3: Δk₂ / k_strict (L = 12 → 18)

- 모든 ε > 0에서 L이 커지면 줄어듭니다. ε ≥ 0.01이면 대략 L^{−1/2}로 줄어 crossover 쪽 거동입니다.
  - g = 0, ε = 0.01: 0.86 → 0.72
  - g = 0.15, ε = 0.0003: 0.991 → 0.985
- ED로는 ε_c에서 Δk, Δe가 0이 되는지 판정할 수 없습니다. REUS/MBAR의 봉우리별 ⟨k⟩, ⟨e⟩와 폭 `sdk_in`, `sdk_act`(`east_mbar`의 `# peaks` 줄, 7장)을 함께 보고 Δk/σ 수준의 봉우리 분리도가 ε와 함께 줄어드는지 검정합니다.

REUS에서 ε = 0.003과 0.01의 χ_KK 최대 위치는 각각 0.0091과 0.0205로, 두 점 기울기는 약 1.63입니다. ED(L = 18)의 같은 두 점 기울기 1.58과 비슷하지만 REUS 위치가 약 0.0022–0.0026 높습니다. 두 REUS 분포는 모두 단일 봉우리이므로 이 값들은 공존점 `s*`가 아니라 pseudocritical crossover 위치이며, 가설 2의 공존선 증거로 사용하지 않습니다.

## 9. 주의사항

- **strict East(ε = 0) TPS는 에르고딕하지 않습니다.** 어느 슬라이스가 빈 격자가 되면 shooting으로 빠져나올 수 없고, 결국 전체가 빈 궤적에 갇힐 수 있습니다. L = 100에서는 평형에서 빈 격자가 나올 확률이 약 10⁻¹⁴이라 사실상 문제없지만, 작은 L에서 ε = 0 TPS 결과는 믿지 마세요.
- **`peak`의 gap은 k = 0 섹터 기준**이고, ε가 작으면 χ_KK 피크가 s 격자 간격보다 훨씬 좁습니다. `peak` 모드는 황금분할로 다시 찾습니다.
- **χ_KK 최댓값이 s 격자의 가장자리에 있을 때**(`edge = 1`)는 s\*를 정의하지 않고 그 행의 s\* 이후 물리량을 모두 `nan`으로 출력합니다(`edge`, `ok` 열은 그대로). 표준 오류로 `note: ... s* reported as nan`이 나옵니다.
  - 현재 `results/peak_T1.dat`에서는 ε = 0.03, g = −0.15 (L = 8–14) 네 행이 해당합니다. 이 점들은 격자 안쪽 s ≈ 0.53에 χ_KK ≈ 0.40의 국소 최댓값이 있지만 L에 따라 커지지 않고, 전체 최댓값은 s = −0.5 쪽(활동적 배경)에 있습니다.
  - 이 판정은 `--s` 범위에 따라 달라지므로 서로 비교할 표는 같은 s 범위(예: `-0.5:0.005:0.8`)로 계산하세요.
- **strict East 연속시간 공존선**은 L에 무관하게 e^{−2s\*} = [(1+2g)² − (1−2c)²]/[4c(1−c)]입니다 (주기경계). 이산시간에는 이런 닫힌 식이 없습니다.
- 연속시간과 이산시간 결과는 시간 단위가 달라 값을 직접 비교하지 않고 경향만 비교합니다.

## 10. 파일 구성

```
Makefile
STATUS.md             작업 현황, 남은 일, 서버 작업 순서
AGENTS.md, CLAUDE.md   서버에서 쓰는 AI 도구용 작업 지침 (CLAUDE.md는 AGENTS.md를 불러옴)
inputs/               예시 입력 (tps, clone, tm_finite, tm_asym, mbar)
scripts/run_ed_plan.sh       ED 계산 일괄 실행
scripts/run_ed_fine.sh       가설 분석용 추가 ED 계산 (ε 세분, L = 16, 18)
scripts/run_ed_analysis.sh   east_edan으로 가설 1–3 분석 → results/analysis/
scripts/run_tps_windows.sh   TPS 창 연쇄 실행 + MBAR (한 조합)
scripts/run_tps_scaling.sh   L, t_obs 스케일링용 TPS + MBAR (여러 조합)
scripts/run_umbrella.sh      K 창 umbrella sampling + MBAR + 공존점
scripts/run_reus.sh          K 창 replica-exchange umbrella sampling + MBAR + 공존점 (체크포인트로 이어서 계산)
scripts/submit_reus.sh       run_reus.sh를 SLURM 작업으로 제출 (노드, 코어, 시간 제한, 자동 requeue)
scripts/slurm_reus.sh        SLURM 작업 스크립트 (submit_reus.sh가 씀)
slurm_jobs/, slurm_logs/     제출할 때 생김: 작업별 매개변수, 제출 목록, 작업 로그
results/              ED 계산 결과, results/analysis/ 가설 분석 결과
src/model, basis, tgen, eig_dense, lanczos, cg, ed, cli_util   연속시간 ED
src/rng       xoshiro256** 난수
src/input     키워드 입력 파일 파서
src/dyn       MCS 동역학, 궤적 저장, shooting
src/stats     블록 jackknife
src/tps       TPS 엔진
src/tm        이산시간 전이행렬 정확해
src/clone     cloning
src/mbar      MBAR (umbrella 창 포함), bootstrap
src/peaks     두 봉우리 분리와 공존점(같은 무게) 탐색
src/reus      replica-exchange umbrella sampling
src/east_*.c  실행 프로그램
```
