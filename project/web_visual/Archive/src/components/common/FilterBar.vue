// FilterBar.vue
<template>
    <div class="filter-bar ui basic segment grid">
        <div class="ui form">
            <div class="inline field">
                <label v-if="label" class="large text">{{ label }}</label>
                <input type="text" v-model="filterText" class="three wide column" @keyup.enter="doFilter" :placeholder="placeholder">
                <button class="ui primary button" @click="doFilter">{{ searchBtnLabel }}</button>
                <button class="ui button" @click="resetFilter">{{ resetBtnLabel }}</button>
            </div>
        </div>
    </div>
</template>

<script>
    export default {
        props: {
            placeholder: "",
            label: {
                type: String,
                default() {
                    return "Search for:"
                }
            },
            searchBtnLabel: {
                type: String,
                default() {
                    return "Go"
                }
            },
            resetBtnLabel: {
                type: String,
                default() {
                    return "Reset"
                }
            }
        },
        data () {
            return {
                filterText: ''
            }
        },
        methods: {
            doFilter () {
                console.log('doFilter:', this.filterText)
                this.$events.fire('filter-set', this.filterText)
            },
            resetFilter () {
                this.filterText = ''
                this.$events.fire('filter-reset')
                console.log('resetFilter')
            }
        }
    }
</script>

<style scoped>
    .inline > button {
        font-size: 1rem;
    }

    .large.text {
       font-size: 1.2rem!important;
    }
</style>