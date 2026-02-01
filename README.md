# OceanSandbox — Oil Spill Response Sim

A visual, interactive Rust simulation of ocean currents and an oil spill response scenario.
Paint oil into the water, watch it advect with the flow, and deploy response tools like booms and skimmers.

Another rust simulation of ocean currents with an oil spill scenario. 
You can paint oil into the water and watch it advect with the flow, and deploy response tools like booms and skimmers to try and contain the spill.

![demo](assets/demo.gif)

## Features
- Real-time 2D current field (wind + swirling flow)
- Oil transport (advection + diffusion + weathering/decay)
- Streamlines + particle tracers for flow visualization
- Response tools:
  - **Booms** (place barriers)
  - **Skimmers** (remove oil locally)
  - **Dispersant** (accelerates breakup/decay in an area)
- Live control panel for wind, diffusion, advection, and tool settings

## Controls
- **Left click**: use current tool (spill oil / boom / skimmer / dispersant / erase)
- **Right click + drag**: pan camera
- **Mouse wheel**: zoom

## Run locally
```bash
cargo run --release
```
