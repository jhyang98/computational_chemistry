# AI 작업 지침 (서버에서 이 프로젝트를 이어서 작업하는 도구용)

## 먼저 읽을 것

1. `STATUS.md`: 지금까지 결과, 남은 일, 서버 작업 순서, 끝난 뒤 확인할 항목
2. `README.md`: 프로그램과 스크립트 사용법 (한국어)

## 사용자 선호 (반드시 지킬 것)

- 코드는 영어만 씁니다. 주석과 출력 메시지도 영어이고 한국어는 넣지 않습니다.
- 새 코드의 기본 언어는 C입니다 (C99, 외부 라이브러리 없음, 기존 Makefile에 추가).
- 실행 파일은 상대 경로로 부릅니다 (`./bin/east_reus`, `bash scripts/run_reus.sh`). 절대 경로는 쓰지 않습니다.
- 사용자에게 주는 실행 안내와 보고는 한국어로 씁니다.
- 정보가 부족하면 추측하지 말고 멈춘 뒤 사용자에게 묻습니다.
- 코드, 문서, 설정, 완료된 결과가 바뀌면 검증 후 GitHub `jhyang98/computational_chemistry`의 `soft_east/` 폴더만 갱신합니다. 원격 저장소의 다른 폴더는 수정하지 않습니다.
- 실행 중 계속 바뀌는 결과는 매 진행 주기마다 올리지 않고, 체크포인트·완료 등 일관된 상태가 되었을 때 갱신합니다.
- 반복적인 진행 상황 같은 자잘한 내용은 문서에 남기지 않습니다. 결과, 중요한 오류와 복구 방법, 확정된 결정처럼 다음 작업에 필요한 참고사항만 `STATUS.md`나 `README.md`에 갱신합니다.
- 무인 감시는 `scripts/monitor_reus_codex.sh <작업번호> <매개변수파일>`을 독립 tmux 세션에서 실행합니다. 셸 감시기는 정상 진행에는 개입하지 않고, 로그 정체·비정상 종료·완료 때만 `codex exec`를 호출합니다.
- 무인 감시가 완료 결과의 문서를 갱신한 뒤 `scripts/sync_github_soft_east.sh`로 원격 저장소의 `soft_east/`만 동기화합니다. 감시기와 `codex exec`가 실제로 시작됐는지 제출 직후 확인합니다.

## 사용자와 확정한 것 (바꾸려면 먼저 사용자에게 물을 것)

- 모델 규약 (`STATUS.md` 1장): 무게 e^{−sK−gE}, 뒤집힘 속도, MCS 규칙, t_obs 정의, E는 t_obs+1개 조각 합
- `kcut auto`: 평활한 P(k)의 두 봉우리 사이 최저점이 경계, 평활은 앞뒤 round(n/(2m))개 구간 평균
- ED `peak` 모드: χ_KK 최댓값이 s 격자 가장자리면 nan
- REUS 시작 단계: 창에 못 들어간 replica는 이웃 창 궤적을 복사해 다시 시작

## SLURM 규칙 (사용자가 정한 조건)

- 계산은 모두 SLURM으로 제출합니다. 로그인 노드에서는 빌드와 짧은 테스트(`make`, `make test`)만 합니다.
  - `make validate`(약 2분)나 다른 계산은 `srun --nodelist=<NODE> --cpus-per-task=4 --time=00:30:00 --export=ALL,OMP_NUM_THREADS=4 make validate OMP=1`처럼 계산노드에서 돌립니다. OpenMP 빌드에서는 할당 코어 수와 같은 `OMP_NUM_THREADS`를 반드시 지정합니다.
- 계산노드 3만 씁니다. `scripts/submit_reus.sh`의 `NODE`가 비어 있으면 `sinfo -N -l`로 노드 이름을 찾고, 사용자에게 맞는지 확인받은 뒤 채웁니다.
- 작업당 `CPUS=4`를 사용하고 동시에 합계 최대 12코어, 즉 REUS 세 작업까지 허용합니다. `srun`으로 다른 계산을 돌릴 때도 같은 시각에 쓰는 코어 합계가 12개를 넘지 않게 합니다.
- 동시에 실행 중인 모든 작업이 사용하는 노드 수의 합은 4개를 넘기지 않습니다. 이 프로젝트의 REUS는 작업마다 node3 한 대만 요청하므로 동시에 최대 세 작업입니다.
- `scripts/submit_reus.sh`는 아직 `JOB_NAME=soft_east`와 singleton을 기본으로 사용합니다. 실행 중인 작업 때문에 스크립트를 바꾸지 말고, 현재 캠페인의 작업 24410·24411은 사용자가 동시 12코어를 허용한 뒤 SLURM dependency만 제거했습니다.
- 한 번 실행은 최대 24시간(`TIME=24:00:00`)입니다. 제출할 때 선택한 SLURM 파티션의 최대 시간이 더 짧으면 그 시간을 자동으로 사용합니다. 안 끝나면 REUS가 체크포인트를 남기고 작업이 자동으로 requeue됩니다(`scripts/slurm_reus.sh`).
- REUS는 항상 `bash scripts/submit_reus.sh [KEY=VALUE …]`로 제출합니다. `scripts/run_reus.sh` 파일을 고쳐서 매개변수를 바꾸지 않습니다.
- 확인 명령
  - `squeue -u $USER -n soft_east`: 대기/실행 중인 작업
  - `tail -n 30 slurm_logs/soft_east_<작업번호>.out`: 진행 상황과 requeue 기록
  - `sacct -j <작업번호> --format=JobID,State,Elapsed,ExitCode,NodeList`
  - 이 서버처럼 SLURM accounting이 꺼져 `sacct`가 실패하면 `scontrol show job <작업번호>`의 `JobState`, `ExitCode`, `Restarts`로 확인
  - `cat slurm_jobs/jobs.txt`: 제출 목록 (작업번호, 결과 폴더, 매개변수 파일)
- 처음 제출한 작업에서 requeue가 실제로 되는지(로그의 `requeue job …` 다음에 같은 작업번호로 다시 시작하는지) 확인하고 사용자에게 알립니다. 클러스터가 사용자 requeue를 막으면 사용자에게 묻습니다.

## 서버 작업 규칙

- 실행 중인 작업이 쓰는 `scripts/`, `bin/`, 결과 폴더는 덮어쓰거나 지우지 않습니다. bash는 실행 중에도 스크립트 파일을 읽으므로, 도중에 파일이 바뀌면 이상하게 동작할 수 있습니다.
- 결과 폴더(`reus/`, `umb/`, `results/`)는 지우지 않습니다. 같은 조건으로 다시 돌릴 때는 `OUT_ROOT`를 바꾸거나 사용자에게 확인합니다.
- 작업 스케줄러(큐)의 종류와 제출 스크립트는 정해 두지 않았습니다. 처음에 사용자에게 확인합니다.
- 실패나 표본 품질 문제를 고치기 위한 REUS 재실행은 조치안과 예상 시간을 기록한 뒤 사용자 동의를 다시 묻지 않고 제출합니다. 끝난 결과 폴더는 덮어쓰지 않고 새 `OUT_ROOT`를 씁니다.
  - 새 (ε, g, L, t_obs) 조합처럼 과학적 조건을 추가하는 계산은 사용자에게 먼저 확인합니다.
  - 분석만 다시 하는 일(`bash scripts/submit_reus.sh MODE=mbar <같은 매개변수>`)은 바로 해도 됩니다.
  - 이미 제출된 작업의 이어서 계산(자동 requeue, `--from`으로 다시 제출)은 동의 없이 해도 됩니다.
- 코드를 고친 뒤에는 항상 다음을 확인하고, 실패하면 결과를 쓰지 않습니다.
  - `make`: 경고 없이 빌드
  - `make test`: 826개 통과 (약 15초)
  - 시뮬레이션·MBAR 코드를 고쳤으면 `make validate`: 89개 통과 (약 2분)
  - 테스트를 추가했으면 위 개수와 `README.md`, `STATUS.md`의 개수를 함께 고칩니다.
- 결과나 동작이 바뀌면 `STATUS.md`(현황)와 `README.md`(사용법)를 한국어로 갱신합니다.

## REUS 모니터링 절차

실행 중 (`slurm_logs/soft_east_<작업번호>.out`)
- `PROGRESS` 주기마다 `cycle …: round trips so far …`가 나옵니다.
- `checkpoint written at cycle …; stopping` → `requeue job …` → `resumed from … at cycle …`: 24시간 제한에 따른 정상 동작입니다.
- `not requeued: already requeued …`: requeue 한도입니다. `bash scripts/submit_reus.sh --from <그 작업의 .env>`로 다시 제출하면 이어집니다.
- `note: window … restarting from the trajectory of window …`: 정상 동작(이웃 궤적 복사)입니다.
- `error: window … not reached …`: 이웃 궤적에서도 못 들어간 것입니다. 작업이 멈춥니다.
  - 조치안: `N_APPROACH` 늘리기, 또는 `INIT`에 예전 궤적 패턴(`umb/…/w%02d.traj`) 지정

끝난 뒤 (`STATUS.md` 5장 표 순서대로)
1. `summary.dat`: round trips, 쌍별 `swap_acc`, `KS`, `n1`·`n2`, 창별 `N_eff`, `from`
2. `mbar.dat` 머리말: MBAR 수렴, 창별 `<K>_raw`와 `<K>_mbar` 일치
3. `mbar_hist.dat`: 공존점 근처 두 봉우리, 활동 봉우리가 `KMAX_ALL` 전에 0으로 감
4. `grep -E "s_equal|^# peaks " mbar.dat`: s\*, 오차, 복제본 수, 봉우리별 k, e

사용자 보고 (한국어, 짧게)
- s\* ± 오차, Δk ± 오차, Δe ± 오차
- 품질 지표 요약: round trips, 가장 나쁜 swap_acc와 KS 쌍
- 문제와 조치안. 오래 걸리는 조치는 동의를 받고 실행합니다.

## 알려진 함정

- 보통 TPS는 L ≥ 50, t_obs ≥ 400에서 이력 현상이 있습니다. 큰 계는 REUS를 씁니다.
- ε = 0(strict East)에서는 빈 궤적에서 빠져나올 수 없습니다. 작은 L의 ε = 0 TPS 결과는 믿지 않습니다.
- 유한 t_obs에서는 비활동 봉우리 궤적도 시간 경계 근처가 활동적이라 k_in이 0보다 꽤 큽니다. Δk, Δe는 t_obs에 따라 비교해야 합니다.
- 분포가 한 봉우리여도 경계로 나누면 dk > 0입니다(정규분포면 약 1.6σ). dk는 sdk_in, sdk_act와 같이 봅니다.
- `east_mbar`의 예전 키워드 `smooth`는 `smooth_res`로 바뀌었습니다.
- 이산시간(MCS)과 연속시간(ED)은 시간 단위가 달라 값은 대략적으로만 비교합니다.
